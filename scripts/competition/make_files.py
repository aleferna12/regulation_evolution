import argparse
import logging
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
    if len(colors) - 1 != len(attrdf):
        raise ValueError("number of non-white pixel colors in 'imgfile' and entries in 'cellfile' "
                         f"table must match (got {len(colors - 1)} and {len(attrdf)} "
                         "respectively)")

    # Reverse colors and indexes so they are ordered by reverse rgb when matching the cells
    colors = colors[::-1]
    lat_size = img.shape[0]
    indexes = np.reshape(indexes, (lat_size, lat_size))
    indexes = len(colors) - 1 - indexes

    # Basically we sample cell_length-interspaced points from image
    group_map = indexes[::cell_length, ::cell_length]
    # Then create a unique sigma for each point where a non-zero value was found in the image
    sigma_lat = np.arange(1, group_map.shape[0]**2 + 1).reshape(group_map.shape)
    sigma_count = np.unique(np.where(group_map, sigma_lat, 0), return_inverse=True)
    sigma_map = sigma_count[1].reshape(sigma_lat.shape)
    # Then tile the pattern to get back to lattice size
    lat = sigma_map.repeat(cell_length, axis=0).repeat(cell_length, axis=1)
    if lat_size % cell_length:
        margin = cell_length - (lat_size % cell_length)
        lat = lat[:-margin, :-margin]

    cells = []
    pbar = Counter(desc="Cells created", total=len(sigma_count[0]) - 1)
    for sigma, group in np.stack([sigma_map, group_map], axis=2).reshape((-1, 2)):
        if sigma != 0:
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


if __name__ == "__main__":
    main()