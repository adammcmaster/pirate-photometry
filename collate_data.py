from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy.time import Time

from tqdm.notebook import tqdm

from pathlib import Path

import itertools
import pickle

from targets import ESCAPED_TARGET_NAMES
import matching
from calibration import calibrate


DATA_PATH = Path("//stem-linux-homes/OSL-Telescope/data/users/Pipeline/")
REPROCESS = False # Set to True to reprocess previous data rather than finding new data

try: 
    with open('data/processed_dates.pickle', 'rb') as processed_dates_file:
        processed_dates = pickle.load(processed_dates_file)
except FileNotFoundError:
    processed_dates = {}
    
if REPROCESS:
    processed_dates = [ d for d, has_data in processed_dates.items() if not has_data ]
    
new_dates = [p for p in DATA_PATH.glob('*/202?_??_??') if p not in processed_dates]
total_dates = len(new_dates)
obs_tables = {}

for date in tqdm(new_dates, desc='New dates'):
    date_has_data = False
    for obs_catalogue in tqdm(list(date.glob('Catalogues/*_anm83_*.cat')), desc=str(date)):
        date_has_data = True
        name = None
        for escaped_target_name, target_name in ESCAPED_TARGET_NAMES.items():
            if escaped_target_name in obs_catalogue.stem:
                name = target_name
                obs_meta = list(itertools.chain.from_iterable(l.split('_') for l in obs_catalogue.stem.split(f'_{escaped_target_name}_')))
                break
        if name is None:
            continue
        obs_meta = {
            'telescope': obs_meta[0],
            'main target': name,
            'band': obs_meta[6][0],
            'exposure': float(obs_meta[6][1:]),
            'timestamp': Time(
                dict(zip(
                    ['year', 'month', 'day', 'hour', 'minute', 'second'],
                    map(int, obs_meta[8:14])
                )),
                format='ymdhms',
            ).jd,
        }

        try:
            Path(f"data/obs_catalogues/").mkdir(parents=True, exist_ok=True)
            table = Table.read(f'data/obs_catalogues/{obs_catalogue.stem}.ecsv')
        except FileNotFoundError:
            table = Table.read(obs_catalogue, format='ascii.sextractor')
            table.write(f'data/obs_catalogues/{obs_catalogue.stem}.ecsv')
        table.rename_column('ALPHA_J2000', 'RA')
        table.rename_column('DELTA_J2000', 'Dec')

        matching.match_targets(table)
        table.rename_column('RA', '_RAJ2000')
        table.rename_column('Dec', '_DEJ2000')
        matched_targets = table[~table['matched target'].mask]
        
        for target_row in matched_targets:
            out_table = calibrate(target_row, obs_meta, table)
 
            for key, val in obs_meta.items():
                out_table[key] = val

            if target_row['matched target'] not in obs_tables:
                try:
                    if REPROCESS:
                        obs_tables[target_row['matched target']] = out_table
                        continue
                    else:
                        obs_tables[target_row['matched target']] = Table.read(f"data/{target_row['matched target']}.ecsv")
                except FileNotFoundError:
                    obs_tables[target_row['matched target']] = out_table
                    continue
            obs_tables[target_row['matched target']] = vstack([obs_tables[target_row['matched target']], out_table])
        
    if not REPROCESS:
        processed_dates[date] = date_has_data

for name, table in obs_tables.items():
    table.write(f"data/{name}.ecsv", overwrite=True)

if not REPROCESS:
    with open('data/processed_dates.pickle', 'wb') as processed_dates_file:
        pickle.dump(processed_dates, processed_dates_file)