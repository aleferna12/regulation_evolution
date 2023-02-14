import argparse
import logging
import pandas as pd
from typing import List
from scripts.analyses.fileio import *

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="make_templates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="build a template file for competition experiments by sampling a single "
                    "cell from each of the input files. Cells are chosen according to how close"
                    "they are to having median food"
    )
    parser.add_argument("cellfiles",
                        nargs='+',
                        help="CSV files containing cell data from simulations")
    parser.add_argument("outfile",
                        help="CSV output template file")
    parser.add_argument("-f",
                        "--food",
                        default=200,
                        type=int,
                        help="the amount of food the will be reassigned to each cell in the "
                             "final template. Use -1 to select the amount of the first "
                             "cell in the final template")
    args = parser.parse_args()
    celldfs = [parse_cell_data(path) for path in args.cellfiles]
    sample_template(celldfs, args.outfile, args.food)


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
    tdf.food = food if food != -1 else tdf.food.iloc[0]
    tdf.time = 0
    tdf.time_since_birth = 0
    tdf.last_meal = 0
    tdf.times_divided = 0
    tdf.dividecounter = 0
    tdf.ancestor = tdf.sigma
    tdf.to_csv(outputfile, index=False)


if __name__ == "__main__":
    main()