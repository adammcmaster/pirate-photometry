from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.coordinates.name_resolve import NameResolveError
from astropy.table import Table, MaskedColumn
from astropy import units as u

import pickle

import targets


try:
    with open("data/target_coords.pickle", "rb") as coords_file:
        target_coords = pickle.load(coords_file)
except (FileNotFoundError, EOFError):
    target_coords = {}

for target_name in targets.TARGETS:
    if target_name not in target_coords:
        try:
            target_coords[target_name] = SkyCoord.from_name(target_name, parse=True)
        except NameResolveError:
            continue

with open("data/target_coords.pickle", "wb") as coords_file:
    pickle.dump(target_coords, coords_file)

target_coords = Table(
    {
        "target_name": list(target_coords.keys()),
        "coords": SkyCoord(list(target_coords.values())),
    }
)


def match_targets(table):
    if "coords" not in table.colnames:
        table["coords"] = SkyCoord.guess_from_table(table)
    nearest_matches, nearest_separations, nearest_distances = match_coordinates_sky(
        table["coords"], target_coords["coords"]
    )
    mask = nearest_separations > (1e-3 * u.deg)
    table["matched target"] = MaskedColumn(
        tuple(target_coords[n]["target_name"] for n in nearest_matches), mask=mask,
    )
    table["matched target separation"] = MaskedColumn(nearest_separations, mask=mask,)

