import argparse
import logging
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

import numpy as np
from enlighten import Counter
from colorir import *
from make_netgraphs import *

config.DEFAULT_COLOR_FORMAT = ColorFormat(sRGB, max_rgb=1, round_to=-1, include_a=True)

logger = logging.getLogger(__name__)
_cs = Palette.load()


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="plot_netgraph",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("datafile",
                        help="either a pickle or a CSV file to read the networks from")
    parser.add_argument("outputdir", help="directory where to save the SVGs")
    parser.add_argument("-p",
                        "--prune-level",
                        help="prune level, the highest the less edges will be shown",
                        default=1,
                        type=int)
    parser.add_argument("-s",
                        "--sigmas",
                        help="list of space-delimited sigmas to plot",
                        default=None,
                        nargs='*')
    args = parser.parse_args()

    if args.datafile[-4:].lower() == ".csv":
        celldf = parse_cell_data(args.datafile)
        netgraphs = make_netgraphs(celldf)
    else:
        netgraphs = read_netgraphs_file(args.datafile)

    if args.sigmas is not None:
        netgraphs = [ng for ng in netgraphs if ng.graph["sigma"]
                     in np.array(args.sigmas, dtype=int)]

    pbar = Counter(total=len(netgraphs), desc="Networks plotted")
    for i, netgraph in enumerate(netgraphs):
        plot_netgraph(netgraph, f"{args.outputdir}/netgraph{i}.svg", args.prune_level)
        pbar.update()

    logger.info("Finished")


def plot_netgraph(netgraph, outfile, prune_level=1):
    # Lowest weight is red, highest is green
    grad = Grad([_cs.red] + [(0, 0, 0, 0)] * prune_level + [_cs.green])
    weights = list(nx.get_edge_attributes(netgraph, "weight").values())
    # To keep 0 at the center we need to have a symmetric distribution
    limit = max(max(weights), -min(weights))
    # Normalize to [0, 1] interval
    edge_colors = [grad.perc((w + limit) / (2 * limit)) for w in weights]
    pos = nx.multipartite_layout(netgraph, subset_key="ntype")

    fig, ax = plt.subplots()
    nx.draw(netgraph, pos, ax=ax, edge_color=edge_colors, with_labels=True)
    fig.savefig(outfile)


if __name__ == "__main__":
    main()