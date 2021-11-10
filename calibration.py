from astropy.coordinates import match_coordinates_sky
from astropy.table import Table, unique
from astropy import units as u
from astroquery.vizier import Vizier

import numpy

from pathlib import Path

MAX_CALIBRATION_ROWS = 10
MIN_CALIBRATION_ROWS = 5
SPLIT_SIZE = 10
CALIBRATION_CATALOGUES = {
    'I/284/out': {
        'B': 'B1mag',
        'R': 'R1mag',
        'I': 'Imag',
    },
    'I/297/out': {
        'B': 'Bmag',
        'V': 'Vmag',
        'R': 'Rmag',
    },
    'VI/135/table15': {
        'B': 'BTmag',
        'V': 'VTmag',
    },
    'II/339/uvotssc1': {
        'U': 'Umag',
        'B': 'Bmag',
        'V': 'Vmag',
    }
}

def get_refstar_table(target_row, obs_meta, obs_table, calibration_catalogue, calibration_catalogue_fields):
    # To do: We need to be smart about caching these when they have fewer than the min rows.
    #        Don't want to repeatedly look up ones which always fail, but don't want to
    #        ignore ones which failed the first time but would succeed later (i.e. with a better
    #        exposure time).
    #        Also when we have multiple cached tables for a given band, it would be good to choose
    #        whichever one has the most rows.
    obs_table = Table(obs_table)
    refstar_table_path = f"data/calibration_tables/{target_row['matched target']}/{calibration_catalogue.replace('/', '_')}_{obs_meta['band']}.ecsv"
    try:
        refstar_table = Table.read(refstar_table_path)
    except FileNotFoundError:
        # Generate a list of reference stars
        obs_table = obs_table[obs_table['FLAGS'] == 0]
        obs_table['flux_diff'] = numpy.abs(obs_table['FLUX_AUTO'] - target_row['FLUX_AUTO'])
        obs_table.sort('flux_diff')
        obs_table = obs_table[:50]
        obs_table['separation'] = obs_table['coords'].separation(target_row['coords'])
        obs_table.sort('separation')

        refstar_table = {
            'coords': [],
            'mag': [],
        }

        total_iterations = int(len(obs_table) / SPLIT_SIZE) + 1
        for i in range(total_iterations):
            sources = obs_table[i * SPLIT_SIZE : (i+1) * SPLIT_SIZE]
            if len(sources) == 0:
                continue

            # Exclude the target star 
            sources = sources[sources['NUMBER'] != target_row['NUMBER']]

            catalogue_results = Vizier.query_region(sources, radius=1e-4*u.deg, catalog=calibration_catalogue)
            if len(catalogue_results) == 0:
                continue

            # Exclude sources with multiple matches
            catalogue_matches = unique(catalogue_results[0], '_q', keep='none')

            # Exclude sources with no mag in the current band
            catalogue_matches = catalogue_matches[~numpy.isnan(catalogue_matches[calibration_catalogue_fields[obs_meta['band']]])]

            # To do: Reject any known variables

            for catalogue_row in catalogue_matches:
                refstar_table['coords'].append(sources[int(catalogue_row['_q']) - 1]['coords'])
                refstar_table['mag'].append(tuple(catalogue_row[[calibration_catalogue_fields[obs_meta['band']]]])[0])

            if len(refstar_table['coords']) >= MAX_CALIBRATION_ROWS:
                break
        refstar_table['coords'] = refstar_table['coords'][:MAX_CALIBRATION_ROWS]
        refstar_table['mag'] = refstar_table['mag'][:MAX_CALIBRATION_ROWS]
        refstar_table = Table(refstar_table)
        if len(refstar_table) < MIN_CALIBRATION_ROWS:
            return None
        refstar_table.write(
            refstar_table_path,
            overwrite=True,
        )
    return refstar_table


def match_refstars(obs_table, refstar_table):
    calibration_table = Table(obs_table)
    nearest_matches, nearest_separations, nearest_distances = match_coordinates_sky(obs_table['coords'], refstar_table['coords'])
    mask = nearest_separations > (1e-3 * u.deg)
    calibration_table['ref mag'] = tuple(refstar_table[n]['mag'] for n in nearest_matches)
    calibration_table['ref star separation'] = nearest_separations
    return calibration_table[~mask]


# To do: How to handle generating the calibration table if the first observation isn't
#        one where this target was the primary target? I.e. the target may be near the
#        edge of the frame.
def calibrate(target_row, obs_meta, obs_table):
    out_table = Table(target_row)

    calibration_succeeded = False

    Path(f"data/calibration_tables/{target_row['matched target']}").mkdir(parents=True, exist_ok=True)
    for calibration_catalogue, calibration_catalogue_fields in CALIBRATION_CATALOGUES.items():
        if obs_meta['band'] not in calibration_catalogue_fields:
            continue
        refstar_table = get_refstar_table(target_row, obs_meta, obs_table, calibration_catalogue, calibration_catalogue_fields)
        if refstar_table is None:
            continue
        
        # Match obs_table to refstar_table
        calibration_table = match_refstars(obs_table, refstar_table)
        if len(calibration_table) < MIN_CALIBRATION_ROWS:
            continue

        # Calculate calibrated mag of target_row

        calibrated_mags = (
            calibration_table['ref mag'] - 2.5 * numpy.log10(target_row['FLUX_AUTO'] / calibration_table['FLUX_AUTO'])
        )

        # To do: Maybe this should be a weighted average using the flux error?
        out_table.add_columns(
            [
                numpy.mean(calibrated_mags),
                # Standard error as per https://www.statology.org/standard-error-of-mean-python/
                numpy.std(calibrated_mags, ddof=1) / numpy.sqrt(numpy.size(calibrated_mags)),
                calibration_catalogue,
            ],
            names=[
                'calibrated_mag',
                'calibrated_mag_err',
                'calibration_catalogue',
            ],
        )
        calibration_succeeded = True
        break

    if not calibration_succeeded:
        out_table.add_columns(
            [
                numpy.nan,
                numpy.nan,
                '',
            ],
            names=[
                'calibrated_mag',
                'calibrated_mag_err',
                'calibration_catalogue',
            ],
        )
    
    return out_table