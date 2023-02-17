# TODO: Update to read CSV files
import logging
import pickle
import networkx as nx
import pandas as pd
import numpy as np
import scripts.fileio as fio
from argparse import ArgumentParser
from pathlib import Path
from enlighten import Counter

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        celldf = fio.parse_cell_data(args.datafile)
        netgraphs = make_netgraphs(celldf)
        with open(args.outputfile, "wb") as file:
            pickle.dump(netgraphs, file)

    parser = ArgumentParser(description="Create network files to be plotted with 'plot_netgraph'")
    parser.add_argument("datafile", help="Input CSV file")
    parser.add_argument("outputfile", help="Output pickle file containing the network graphs")
    parser.set_defaults(run=run)
    return parser


def read_netgraphs_file(filepath):
    with open(filepath, "rb") as file:
        return pickle.load(file)


def make_netgraphs(celldf):
    logger.info("Creating graphs from the cells' GRNs")

    netgraphs = []
    pbar = Counter(total=len(celldf.index), desc="Networks parsed")
    for index in celldf.index:
        netgraphs.append(make_netgraph(celldf.loc[index]))
        pbar.update()

    return netgraphs


def make_netgraph(cellss: pd.Series):
    """Make directed network graph from cell series."""
    cellss = cellss.copy()
    for key in fio._gene_attrs:
        cellss[key] = np.fromstring(cellss[key], sep=' ')

    netgraph = nx.DiGraph(sigma=cellss["sigma"])

    # ntype is used to plot and must be numerical, but is the same as type
    for i in range(cellss["innr"]):
        netgraph.add_node(i, type="in", ntype=0, scale=cellss["in_scale_list"][i])

    for i in range(cellss["regnr"]):
        reg_id = cellss["innr"] + i
        netgraph.add_node(reg_id, type="reg", ntype=1, threshold=cellss["reg_threshold_list"][i])
        for j in range(cellss["innr"]):
            w_index = i * cellss["innr"] + j
            netgraph.add_edge(j, reg_id, weight=cellss["reg_w_innode_list"][w_index])
        for j in range(cellss["regnr"]):
            reg_id2 = cellss["innr"] + j
            w_index = i * cellss["regnr"] + j
            netgraph.add_edge(reg_id2, reg_id, weight=cellss["reg_w_regnode_list"][w_index])

    for i in range(cellss["outnr"]):
        out_id = cellss["innr"] + cellss["regnr"] + i
        netgraph.add_node(out_id, type="out", ntype=2, threshold=cellss["out_threshold_list"][i])
        for j in range(cellss["regnr"]):
            reg_id = cellss["innr"] + j
            w_index = i * cellss["regnr"] + j
            netgraph.add_edge(reg_id, out_id, weight=cellss["out_w_regnode_list"][w_index])

    return netgraph
