from astropy.coordinates import match_coordinates_sky
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import Table, unique, vstack
from astropy import units as u
from astroquery.vizier import Vizier

import numpy

from tqdm import tqdm

import constants
from matching import target_coords
from targets import (
    get_target_observations,
    read_source_catalogue,
    TARGETS,
    write_target_observations,
)

MAX_CALIBRATION_ROWS = 10
MIN_CALIBRATION_ROWS = 5
SPLIT_SIZE = 10
CALIBRATION_CATALOGUES = {
    "I/284/out": {"B": "B1mag", "R": "R1mag", "I": "Imag"},
    "I/297/out": {"B": "Bmag", "V": "Vmag", "R": "Rmag"},
    "VI/135/table15": {"B": "BTmag", "V": "VTmag"},
    "II/339/uvotssc1": {"U": "Umag", "B": "Bmag", "V": "Vmag"},
}
CALIBRATION_CATALOGUES_PER_BAND = {}
for cat, bands in CALIBRATION_CATALOGUES.items():
    for band in bands:
        CALIBRATION_CATALOGUES_PER_BAND.setdefault(band, set()).add(cat)

CATALOGUE_MATCHING_RADIUS = 1e-4 * u.deg


def calibrate_data():
    for target in tqdm(TARGETS.keys(), desc="Calibrating observations"):
        refstars = get_refstars_for_target(target)
        if refstars is None:
            continue
        target_obs = get_target_observations(target)

        mags = []
        errors = []
        for row in target_obs:
            band_refstars = refstars[refstars["band"] == row["band"]]
            if len(band_refstars) == 0:
                mags.append(None)
                errors.append(None)
                continue

            # Skip this row if the mag and error have already been calculated
            try:
                mag = row["calibrated magnitude"]
                err = row["calibrated magnitude err"]
                if (
                    mag is not None and err is not None
                    and not numpy.ma.is_masked(mag) 
                    and not numpy.ma.is_masked(err)
                ):
                    mags.append(mag)
                    errors.append(err)
                    continue
            except KeyError:
                pass

            source_catalogue = read_source_catalogue(row["observation catalogue"])
            obs_refstars = match_refstars(
                source_catalogue, refstars[refstars["band"] == row["band"]]
            )
            mag, err = calibrate(row, obs_refstars)
            mags.append(mag)
            errors.append(err)

        target_obs["calibrated magnitude"] = mags
        target_obs["calibrated magnitude err"] = errors
        write_target_observations(target, target_obs)


def get_refstars_for_target(target):
    """
    Generate a catalogue of reference stars for a given target.

    - target: The target name.

    Returns: An astropy table of ref star coordinates and magnitudes.
    """
    refstars_path = constants.REFERENCE_STARS_PATH / f"{target}.ecsv"
    if refstars_path.exists():
        return Table.read(refstars_path)

    out_refstars = []

    # Select up to 10 unflagged observations of the target in each band at random.
    band_catalogues = get_random_catalogues_for_bands(target)
    if band_catalogues is None:
        return

    for band, (catalogues, target_rows) in band_catalogues.items():
        # Cut original observation catalogues on angular separation and absolute flux difference.
        ref_stars = cut_catalogues_flux(vstack(target_rows)["FLUX_AUTO"], catalogues)
        ref_stars = cut_catalogues_separation(
            target_coords[target_coords["target_name"] == target]["coords"], ref_stars
        )

        # Concatenate the original observation catalogues into one table.
        ref_stars = vstack(ref_stars)

        # Cross-reference the candidates against each catalogue, keeping whichever has
        # the most matches
        ref_stars = get_most_catalogue_matches(ref_stars, band)
        if ref_stars is None:
            continue

        # Remove known variable stars
        ref_stars = remove_known_variables(ref_stars)

        # To do: Remove stars which show variability in our observations

        ref_stars["band"] = band
        out_refstars.append(ref_stars)

    if len(out_refstars) > 0:
        out_refstars = vstack(out_refstars)
        out_refstars["RA"] = out_refstars["_RAJ2000"]
        out_refstars["Dec"] = out_refstars["_DEJ2000"]
        out_refstars["coords"] = SkyCoord.guess_from_table(out_refstars)
        out_refstars.write(refstars_path, overwrite=True)
        return out_refstars


def get_random_catalogues_for_bands(target, catalogues_per_band=10):
    """
    Read source catalogues for a random set of unflagged observations of the target.

    - target: The target name.

    Returns: A dict of lists of tuples per band with full catalogue tables and 
             individual target observations:
             {'R': ([Table(), Table(), ...], [Row(), Row(), ...]), ...}
    """
    target_obs = get_target_observations(target)
    if target_obs is None:
        return
    target_obs = target_obs[
        target_obs["FLAGS"] == 0
    ]  # To do: Should this only care about saturation?

    out = {}
    for band_obs in target_obs.group_by("band").groups:
        num_band_obs = len(band_obs)
        if num_band_obs < catalogues_per_band:
            continue
        # Select a random subset of observations
        selected_obs = [
            row
            for row in band_obs[
                numpy.random.choice(num_band_obs, catalogues_per_band, replace=False)
            ]
        ]
        # Load the full source tables for each selected observation
        selected_tables = [
            read_source_catalogue(row["observation catalogue"]) for row in selected_obs
        ]
        # Remove the target's observation from the full source table
        selected_tables = [
            table[table["NUMBER"] != obs_row["NUMBER"]]
            for table, obs_row in zip(selected_tables, selected_obs)
        ]
        out[selected_obs[0]["band"]] = (selected_tables, selected_obs)

    return out


def cut_catalogues_flux(fluxes, catalogues, keep=50):
    """
    Selects sources by smallest flux difference in a list of Tables.

    - fluxes: A list of absolute flux differences
    - catalogues: [Table(), Table(), ...]
    - keep: the number of sources to keep in each table

    Returns: a list of filtered Tables
    """
    out_tables = []
    for flux, table in zip(fluxes, catalogues):
        table["flux diff"] = numpy.abs(flux - table["FLUX_AUTO"])
        table.sort("flux diff")
        out_tables.append(table[:keep])
    return out_tables


def cut_catalogues_separation(coords, catalogues, keep=25):
    """
    Selects sources by closest angular separation in a dict of lists of Tables.

    - coords: The coordinates of the target
    - catalogues: [Table(), Table(), ...]
    - keep: the number of sources to keep in each table

    Returns: a list of filtered Tables.
    """
    out_tables = []
    for table in catalogues:
        table["separation"] = SkyCoord.guess_from_table(table).separation(coords)
        table.sort("separation")
        out_tables.append(table[:keep])
    return out_tables


def get_most_catalogue_matches(stars, band):
    most_matches = 0
    best_matches = None
    for catalogue in CALIBRATION_CATALOGUES_PER_BAND[band]:
        catalogue_results = Vizier.query_region(
            stars, radius=CATALOGUE_MATCHING_RADIUS, catalog=catalogue
        )
        if len(catalogue_results) == 0:
            continue
        else:
            catalogue_results = catalogue_results[0]

        # Remove multiple matches
        catalogue_results = unique(catalogue_results, "_q", keep="none")

        # Remove results with no magnitude
        catalogue_results = catalogue_results[
            ~numpy.isnan(catalogue_results[CALIBRATION_CATALOGUES[catalogue][band]])
        ]

        results_len = len(catalogue_results)
        if results_len > most_matches:
            most_matches = results_len
            best_matches = catalogue_results

    if best_matches is None:
        return

    unique_matches = {
        "_RAJ2000": [],
        "_DEJ2000": [],
        "magnitude": [],
        "num_matches": [],
    }
    for grouped_matches in best_matches.group_by(["RAJ2000", "DEJ2000"]).groups:
        unique_matches["num_matches"].append(len(grouped_matches))
        unique_matches["_RAJ2000"].append(grouped_matches[0]["RAJ2000"])
        unique_matches["_DEJ2000"].append(grouped_matches[0]["DEJ2000"])
        unique_matches["magnitude"].append(
            grouped_matches[0][CALIBRATION_CATALOGUES[best_matches.meta["name"]][band]]
        )

    return Table(unique_matches, units={"_RAJ2000": u.deg, "_DEJ2000": u.deg})


def remove_known_variables(stars):
    vsx_matches = Vizier.query_region(
        stars, radius=CATALOGUE_MATCHING_RADIUS, catalog="B/vsx/vsx"
    )
    if len(vsx_matches) == 0:
        return stars
    stars.remove_rows(vsx_matches[0]["_q"] - 1)
    return stars


def match_refstars(obs_table, refstar_table):
    obs_table = Table(obs_table)
    nearest_matches, nearest_separations, _ = match_coordinates_sky(
        obs_table["coords"], refstar_table["coords"]
    )
    mask = nearest_separations <= (1e-3 * u.deg)
    obs_table["reference magnitude"] = tuple(
        refstar_table[n]["magnitude"] for n in nearest_matches
    )
    obs_table["reference star separation"] = nearest_separations
    return obs_table[mask]


def calibrate(row, refstars):
    """
    Takes a row and a table of reference stars and returns the calibrated magnitude

    - row: A Row() from a target observation file.
    - refstars: The target's refstars table.

    Returns: A tuple with calibrated magnitude and calibrated magnitude error
    """
    calibrated_mags = refstars["reference magnitude"] - 2.5 * numpy.log10(
        row["FLUX_AUTO"] / refstars["FLUX_AUTO"]
    )

    calibrated_magnitude = numpy.mean(calibrated_mags)
    calibrated_magnitude_err = numpy.std(calibrated_mags, ddof=1) / numpy.sqrt(
        numpy.size(calibrated_mags)
    )

    return calibrated_magnitude, calibrated_magnitude_err
