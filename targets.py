from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import Table

from datetime import date
import re

import constants

ALL_FILTERS = {"B", "V", "R", "I"}

ALL_TELESCOPE_FILTERS = {
    "PIRATE": ALL_FILTERS,
    "COAST": ALL_FILTERS - {"I"},
}

TARGETS = {
    "HD38451": {"end": date(2022, 1, 1)},
    "27Cephei": {"end": date(2021, 11, 9)},
    "1SWASPJ002552.75+454445.3": {
        "expected_periods": (24671.86133,),
        "uses": "",
        "end": date(2021, 12, 15),
    },
    "1SWASPJ002552.88+454558.7": {
        "expected_periods": (148009.10938,),
        "uses": "",
        "end": date(2021, 12, 15),
    },
    "1SWASPJ184559.85+471858.4": {
        "expected_periods": (37009.83594,),
        "uses": "",
        "end": date(2021, 12, 15),
    },
    "1SWASPJ002551.12+454523.8": {
        "expected_periods": (74014.22656,),
        "uses": "",
        "end": date(2021, 12, 15),
    },
    "1SWASPJ054938.69+215431.1": {
        "expected_periods": (20139.55078,),
        "uses": "",
        "end": date(2022, 1, 30),
    },
}

ESCAPED_TARGET_NAMES = {re.sub(r"[.+-]", "_", name): name for name in TARGETS.keys()}


def target_observations_path(target):
    return constants.TARGET_OBSERVATIONS_PATH / f"{target}.ecsv"


def get_target_observations(target, default=None):
    obs_path = target_observations_path(target)
    if obs_path.exists():
        return Table.read(constants.TARGET_OBSERVATIONS_PATH / f"{target}.ecsv")
    else:
        return default


def write_target_observations(target, obs):
    obs.write(target_observations_path(target), overwrite=True)


def read_source_catalogue(path):
    table = Table.read(path, format="ascii.sextractor")
    table.rename_column("ALPHA_J2000", "RA")
    table.rename_column("DELTA_J2000", "Dec")
    table["_RAJ2000"] = table["RA"]
    table["_DEJ2000"] = table["Dec"]
    table["coords"] = SkyCoord.guess_from_table(table)
    return table
