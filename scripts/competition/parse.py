import argparse
import logging
from matplotlib.pyplot import imread

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="parse")
    parser.add_argument(
        "filepath",
        help="image file containing a drawing of the initial setting of the simulation"
    )
