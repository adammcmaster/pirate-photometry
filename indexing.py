from astropy.table import Table, unique, vstack
from astropy.time import Time
from tqdm import tqdm
from pathlib import Path

import itertools
import pickle

from matching import match_targets
from targets import (
    ESCAPED_TARGET_NAMES,
    get_target_observations,
    write_target_observations,
    read_source_catalogue,
)

import constants


# A dict of the dates which have been processed, with boolean values
# indicating whether each one contained any data
INDEXED_DATES_PATH = constants.DATA_PATH / "indexed_dates.pickle"

CATALOGUES_GLOB = "Catalogues/*_anm83_*.cat"


def index_data(reprocess=False):
    """
    Finds data in the osl-telescope Pipeline directory and generates or updates tables with
    observations of each target.

    - reprocess: Whether to skip finding dates with new data and just re-process the previous ones.
    """
    if INDEXED_DATES_PATH.exists():
        with INDEXED_DATES_PATH.open("rb") as processed_dates_file:
            processed_dates = pickle.load(processed_dates_file)
    else:
        processed_dates = {}

    if reprocess:
        processed_dates = [d for d, has_data in processed_dates.items() if not has_data]

    new_dates = [
        p
        for p in constants.OSL_DATA_PATH.glob("*/202?_??_??")
        if p not in processed_dates
    ]
    obs_tables = {}

    if reprocess:
        desc = "Reprocessing indexed dates"
    else:
        desc = "Indexing new dates"

    for date in tqdm(new_dates, desc=desc):
        date_has_data = False
        for obs_catalogue in date.glob(CATALOGUES_GLOB):
            date_has_data = True

            name = None
            for escaped_target_name, target_name in ESCAPED_TARGET_NAMES.items():
                if escaped_target_name in obs_catalogue.stem:
                    name = target_name
                    obs_meta = list(
                        itertools.chain.from_iterable(
                            l.split("_")
                            for l in obs_catalogue.stem.split(
                                f"_{escaped_target_name}_"
                            )
                        )
                    )
                    break
            if name is None:
                continue
            obs_meta = {
                "telescope": obs_meta[0],
                "main target": name,
                "band": obs_meta[6][0],
                "exposure": float(obs_meta[6][1:]),
                "timestamp": Time(
                    dict(
                        zip(
                            ["year", "month", "day", "hour", "minute", "second"],
                            map(int, obs_meta[8:14]),
                        )
                    ),
                    format="ymdhms",
                ).jd,
            }

            table = read_source_catalogue(obs_catalogue)
            match_targets(table)
            matched_targets = table[~table["matched target"].mask]

            for target_row in matched_targets:
                out_table = Table(target_row)

                for key, val in obs_meta.items():
                    out_table[key] = val

                out_table["observation catalogue"] = str(obs_catalogue)

                if target_row["matched target"] in obs_tables:
                    obs_tables[target_row["matched target"]] = vstack(
                        [obs_tables[target_row["matched target"]], out_table]
                    )
                else:
                    if reprocess:
                        obs_tables[target_row["matched target"]] = out_table
                    else:
                        obs_tables[
                            target_row["matched target"]
                        ] = get_target_observations(
                            target_row["matched target"], default=out_table
                        )

        if not reprocess:
            processed_dates[date] = date_has_data

    if not reprocess:
        with INDEXED_DATES_PATH.open("wb") as processed_dates_file:
            pickle.dump(processed_dates, processed_dates_file)

    for name, table in obs_tables.items():
        write_target_observations(name, table)
