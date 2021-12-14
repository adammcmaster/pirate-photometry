from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

import csv
import datetime

from targets import TARGETS, ALL_TELESCOPE_FILTERS, ALL_FILTERS, get_target_observations

import constants


MAX_FLAG_THRESHOLD = 0.25
MIN_EXPOSURES = 5
DEFAULT_EXPOSURES = list(range(1, 5)) + list(range(5, 60, 5)) + list(range(60, 660, 60))
DATE_FORMAT = "%Y-%m-%dT09:00:00"

PRIORITY = 1000
BINNING = 1
REPEATS = 1
MAX_USES = 3
ALT = 20
IS_TIMED = "FALSE"
START_TIME = None
END_TIME = datetime.date.today() + datetime.timedelta(days=21)


def generate_schedules():
    for telescope, _ in ALL_TELESCOPE_FILTERS.items():
        with (constants.DATA_PATH / f"schedule_{telescope}.csv").open("w") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                (
                    "Name",
                    "RA",
                    "Dec",
                    "is timed",
                    "start time",
                    "end time",
                    "Priority",
                    "Atoms",
                    "Binning",
                    "repeats",
                    "Max uses",
                    "Alt limit",
                )
            )
            for target, config in TARGETS.items():
                filters = config.get("filters", ALL_FILTERS)
                filters = filters & ALL_TELESCOPE_FILTERS[telescope]
                _filters = {}
                for f in filters:
                    _filters[f] = []
                    for e in DEFAULT_EXPOSURES:
                        _filters[f].append(float(e))
                filters = _filters

                photometry = get_target_observations(target)
                if photometry is not None:
                    photometry = photometry[photometry["telescope"] == telescope]
                    photometry["flagged"] = (
                        photometry["FLAGS"] & 4
                    ) != 0  # In this case we only care about saturation
                    for band, group in photometry.to_pandas().groupby("band"):
                        if band not in filters:
                            continue
                        pt = group[["exposure", "flagged"]].pivot_table(
                            columns="flagged",
                            index=["exposure"],
                            aggfunc=lambda x: len(x),
                            fill_value=0,
                        )
                        if True not in pt:
                            pt[True] = 0
                        if False not in pt:
                            pt[False] = 0
                        pt["total"] = pt[True] + pt[False]
                        pt["frac"] = pt[True] / pt["total"]

                        flagged_exposures = pt[
                            (pt["total"] >= MIN_EXPOSURES)
                            & (pt["frac"] >= MAX_FLAG_THRESHOLD)
                        ].index
                        if len(flagged_exposures) > 0:
                            min_flagged_exposure = flagged_exposures.min()
                            filters[band] = [
                                e for e in filters[band] if e < min_flagged_exposure
                            ]

                        successful_exposures = pt[
                            (pt["total"] >= MIN_EXPOSURES)
                            & (pt["frac"] < MAX_FLAG_THRESHOLD)
                        ].index
                        if len(successful_exposures) > 0:
                            max_successful_exposure = successful_exposures.max()
                            filters[band] = [
                                e for e in filters[band] if e >= max_successful_exposure
                            ]

                coords = SkyCoord.from_name(target, parse=True)
                atoms = []
                for f, es in filters.items():
                    if len(es) == 0:
                        continue
                    for e in {es[0], es[int(len(es) / 2)], es[-1]}:
                        # Take the unique lowest, middle, and highest exposure times
                        atoms.append(f"{f}:{e}")
                atoms = ";".join(atoms)
                if not atoms:
                    continue

                start = config.get("start", START_TIME)
                if start is None:
                    start = ""
                else:
                    start = start.strftime(DATE_FORMAT)

                end = config.get("end", END_TIME)
                if end is None:
                    end = ""
                else:
                    if end < datetime.date.today():
                        continue
                    end = end.strftime(DATE_FORMAT)

                csv_writer.writerow(
                    [
                        target,
                        coords.ra.to_string(alwayssign=True, sep=":", unit=u.hour),
                        coords.dec.to_string(alwayssign=True, sep=":"),
                        IS_TIMED,
                        start,
                        end,
                        config.get("priority", PRIORITY),
                        atoms,
                        BINNING,
                        config.get("repeats", REPEATS),
                        config.get("uses", MAX_USES),
                        config.get("alt", ALT),
                    ]
                )

