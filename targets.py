from datetime import date
import re

ALL_FILTERS = {"U", "B", "V", "R", "I"}

ALL_TELESCOPE_FILTERS = {
    "PIRATE": ALL_FILTERS,
    "COAST": ALL_FILTERS - {"U", "I"},
}

TARGETS = {
    "HD38451": {"end": date(2022, 1, 1),},
    "27Cephei": {"end": date(2021, 11, 9),},
    "1SWASPJ002552.75+454445.3": {"expected_periods": (24671.86133,), "uses": "",},
    "1SWASPJ002552.88+454558.7": {"expected_periods": (148009.10938,), "uses": "",},
    "1SWASPJ184559.85+471858.4": {"expected_periods": (37009.83594,), "uses": "",},
    "1SWASPJ002551.12+454523.8": {"expected_periods": (74014.22656,), "uses": "",},
}

ESCAPED_TARGET_NAMES = {re.sub(r"[.+-]", "_", name): name for name in TARGETS.keys()}
