from astropy.table import Table

from pathlib import Path

NUM_CALIBRATION_ROWS = 10
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

# To do: How to handle generating the calibration table if the first observation isn't
#        one where this target was the primary target? I.e. the target may be near the
#        edge of the frame.
def calibrate(target_row, obs_meta, table):
    out_table = Table(target_row)

    calibration_succeeded = False

    Path(f"data/calibration_tables/{target_row['matched target']}").mkdir(parents=True, exist_ok=True)
    for calibration_catalogue, calibration_catalogue_fields in CALIBRATION_CATALOGUES.items():
        if obs_meta['band'] in calibration_catalogue_fields:
            calibration_catalogue_path = f"data/calibration_tables/{target_row['matched target']}/{calibration_catalogue.replace('/', '_')}_{obs_meta['band']}.ecsv"
            try:
                calibration_table = Table.read(calibration_catalogue_path)
            except FileNotFoundError:
                table = table[table['FLAGS'] == 0]
                table['flux_diff'] = numpy.abs(table['FLUX_AUTO'] - target_row['FLUX_AUTO'])
                table.sort('flux_diff')
                table = table[:50]
                table['separation'] = table['coords'].separation(target_row['coords'])
                table.sort('separation')

                calibration_rows = []

                total_iterations = int(len(table) / SPLIT_SIZE) + 1
                for i in range(total_iterations):
                    sources = table[i * SPLIT_SIZE : (i+1) * SPLIT_SIZE]
                    if len(sources) == 0:
                        continue
                    sources = vstack([ row for row in sources if row['NUMBER'] != target_row['NUMBER'] ])
                    catalogue_results = Vizier.query_region(sources, radius=1e-4*u.deg, catalog=calibration_catalogue)
                    if len(catalogue_results) == 0:
                        continue

                    # Reject any sources with multiple matches
                    catalogue_matches = unique(catalogue_results[0], '_q', keep='none')

                    # To do: Reject any known variables

                    for catalogue_row in catalogue_matches[~numpy.isnan(catalogue_matches[calibration_catalogue_fields[obs_meta['band']]])]:
                        calibration_rows.append(hstack([
                            sources[int(catalogue_row['_q']) - 1],
                            catalogue_row[[calibration_catalogue_fields[obs_meta['band']]]],
                        ]))

                    if len(calibration_rows) >= NUM_CALIBRATION_ROWS:
                        break
                if len(calibration_rows) < NUM_CALIBRATION_ROWS:
                    continue
                calibration_table = vstack(calibration_rows[:NUM_CALIBRATION_ROWS], metadata_conflicts='silent')
                calibration_table.write(
                    calibration_catalogue_path,
                    overwrite=True,
                )

            calibrated_mags = (
                calibration_table[calibration_catalogue_fields[obs_meta['band']]] 
                - 2.5 * numpy.log10(target_row['FLUX_AUTO'] / calibration_table['FLUX_AUTO'])
            )

            # To do: Maybe this should be a weighted average using the flux error?
            out_table.add_columns(
                [
                    numpy.mean(calibrated_mags),
                    # Standard error as per https://www.statology.org/standard-error-of-mean-python/
                    numpy.std(calibrated_mags, ddof=1) / numpy.sqrt(numpy.size(calibrated_mags)),
                    calibration_catalogue,
                    obs_meta['main target'],
                ],
                names=[
                    'calibrated_mag',
                    'calibrated_mag_err',
                    'calibration_catalogue',
                    'main_target',
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
                obs_meta['main target'],
            ],
            names=[
                'calibrated_mag',
                'calibrated_mag_err',
                'calibration_catalogue',
                'main_target',
            ],
        )
    
    return out_table