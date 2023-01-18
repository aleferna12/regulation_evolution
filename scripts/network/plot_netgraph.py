import argparse
import logging
import networkx as nx
import matplotlib.pyplot as plt
from colorir import *
from make_netgraphs import *

config.DEFAULT_COLOR_FORMAT = ColorFormat(sRGB, max_rgb=1, round_to=-1, include_a=True)

logger = logging.getLogger(__name__)
_cs = Palette.load()


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="plot_netgraph")

    datafile = "/home/aleferna/CProjects/Projects/regulation_evolution/runs/debug/netgraphs.pickle"
    netgraphs = read_netgraphs_file(datafile)
    netgraph = netgraphs[0]

    # Red is -1, green is 1
    grad = Grad([_cs.red, (0, 0, 0, 0), _cs.green])
    weights = list(nx.get_edge_attributes(netgraph, "weight").values())
    print(weights)
    min_w, max_w = min(weights), max(weights)
    # Normalize to [0, 1] interval
    edge_colors = [grad.perc((w - min_w) / (max_w - min_w)) for w in weights]
    pos = nx.multipartite_layout(netgraph)
    nx.draw(netgraph, pos, edge_color=edge_colors, with_labels=True)
    plt.show()


if __name__ == "__main__":
    main()