import argparse
import logging
from typing import List

import numpy as np
import pandas as pd
from matplotlib.pyplot import imread
from enlighten import Counter
from scripts.analyses.fileio import parse_cell_data

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="make_files",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    input_ = parser.add_argument_group("input")
    input_.add_argument(
        "cellfile",
        help="CSV file containing the template cells (number of rows must match number "
             "of colors in the image). The attributes of the template will be copied to all cells"
             "in the same color group (besides sigma, which will be reassigned; and color, which"
             "will default to the group index on the table)"
    )
    input_.add_argument(
        "imgfile",
        help="PNG file containing a drawing of the initial setting of the simulation."
             "Colors represent entries on the 'cellfile' table (ordered in reverse by RGB value)."
             "The background must be white and will be ignored. Be careful not to have any fading "
             "colors, as these will be interpreted as additional cell templates. The image "
             "dimensions should match those of the simulation lattice"
    )
    output = parser.add_argument_group("output")
    output.add_argument("outcellfile", help="CSV output file containing cell data")
    output.add_argument("latticefile", help="CSV output lattice file")
    parser.add_argument(
        "-l",
        "--cell_length",
        help="diameter of the initialized cells",
        default=7,
        type=int
    )
    args = parser.parse_args()

    celldf, latdf = make_competition(args.imgfile, args.cellfile, args.cell_length)
    logger.info("Writing output competition files")
    celldf.to_csv(args.outcellfile, index=False)
    latdf.to_csv(args.latticefile, header=False, index=False)


def make_competition(imgfile, cellfile, cell_length=7):
    logger.info("Making competition file from image and templates")

    if imgfile[-3:].lower() != "png":
        raise ValueError("input image must be a PNG file")
    img = imread(imgfile)
    if img.shape[0] != img.shape[1]:
        raise ValueError("image width and height must match")

    attrdf = parse_cell_data(cellfile)
    attrdf = attrdf.reset_index(drop=True)
    attrdf["time"] = 0

    colors, indexes = np.unique(np.reshape(img, (-1, img.shape[2])), axis=0, return_inverse=True)
    indexes = np.reshape(indexes, (img.shape[0], img.shape[0]))
    # Reverse colors and indexes so they are ordered by reverse rgb when matching the cells
    colors = colors[::-1]
    indexes = len(colors) - 1 - indexes

    if len(colors) - 1 != len(attrdf):
        raise ValueError("number of non-white pixel colors in 'imgfile' and entries in 'cellfile' "
                         f"table must match (got {len(colors - 1)} and {len(attrdf)} "
                         "respectively)")

    sq_lat = make_square_lat(img.shape[0], cell_length)
    # Mask sq_lat with indexes
    lat = np.where(indexes != 0, sq_lat, 0)
    cells = []
    pbar = Counter(total=len(colors) - 1, desc="Groups processed")
    for group in range(1, len(colors)):
        group_cells = lat[indexes == group]
        gpx_count = dict(zip(*np.unique(group_cells, return_counts=True)))
        for sigma, count in gpx_count.items():
            # Only sample whole cells
            if count < cell_length**2:
                continue
            cell_attrs = attrdf.iloc[group - 1].copy()
            cell_attrs["sigma"] = sigma
            cell_attrs["group"] = group - 1
            cells.append(cell_attrs)
        pbar.update()
    celldf = pd.DataFrame(cells).set_index("sigma", drop=False).sort_index()
    celldf.index.name = "sigma_i"

    # Remove cells that are not whole and borders
    trimmed_lat = np.where(np.isin(lat, celldf["sigma"]), lat, 0)[1:-1, 1:-1]
    latdf = pd.DataFrame(trimmed_lat)

    return celldf, latdf


def make_square_lat(lat_size, cell_length):
    # Making sure cells * cells_in_lat = lat_size
    cells_in_lat = lat_size // cell_length + 1
    cells = np.arange(1, cells_in_lat**2 + 1).reshape((cells_in_lat, cells_in_lat))
    lat = np.kron(cells, np.ones((cell_length, cell_length)))
    return lat[:lat_size, :lat_size].astype(int)


def sample_template(celldfs: List[pd.DataFrame], outputfile, food=200):
    """Tries to sample a good cell from cell data frames and resets some attributes."""
    templates = []
    for i, celldf in enumerate(celldfs, 1):
        celldf = celldf[celldf.tau == 1]
        # Sort by closest to median food
        med_food = celldf.food.median()
        template = celldf.iloc[(celldf.food - med_food).abs().argsort()[0]]
        template.sigma = i
        templates.append(template)
    tdf = pd.DataFrame(templates)
    tdf.food = food
    tdf.time = 0
    tdf.time_since_birth = 0
    tdf.last_meal = 0
    tdf.times_divided = 0
    tdf.dividecounter = 0
    tdf.ancestor = tdf.sigma
    tdf.to_csv(outputfile, index=False)


if __name__ == "__main__":
    main()