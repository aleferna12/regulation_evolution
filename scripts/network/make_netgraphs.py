# TODO: Update to read CSV files
import logging
import pickle
import networkx as nx
from argparse import ArgumentParser
from pathlib import Path
from enlighten import Counter
from parse import *

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = ArgumentParser(prog="make_netgraphs")
    parser.add_argument("datafile", help="input CSV file")
    parser.add_argument("outputfile", help="output pickle file containing the network graphs")
    args = parser.parse_args()

    celldf = parse_cell_data(args.datafile)
    netgraphs = make_netgraphs(celldf)

    with open(args.outputfile, "wb") as file:
        pickle.dump(netgraphs, file)


def read_netgraphs_file(filepath):
    with open(filepath, "rb") as file:
        return pickle.load(file)


def make_netgraphs(celldf):
    logger.info("Creating graphs from the cells' GRNs")

    for key in ["in_scale_list",
                "reg_threshold_list",
                "reg_w_innode_list",
                "reg_w_regnode_list",
                "out_threshold_list",
                "out_w_regnode_list"]:
        if isinstance(celldf[key].iat[0], str):
            celldf[key] = celldf[key].apply(lambda s: np.fromstring(s, sep=' '))

    netgraphs = []
    pbar = Counter(total=len(celldf.index), desc="Networks parsed")
    for index in celldf.index:
        netgraphs.append(make_netgraph(celldf.loc[index]))
        pbar.update()

    return netgraphs


def make_netgraph(cellss):
    """Make directed network graph from cell series.
    The string attributes need to have been already parsed into lists of floats.
    """
    netgraph = nx.DiGraph()

    for i in range(cellss["innr"]):
        netgraph.add_node(i, type="in", subset=0, scale=cellss["in_scale_list"][i])

    for i in range(cellss["regnr"]):
        reg_id = cellss["innr"] + i
        netgraph.add_node(reg_id, type="reg", subset=1, threshold=cellss["reg_threshold_list"][i])
        for j in range(cellss["innr"]):
            w_index = i * cellss["innr"] + j
            netgraph.add_edge(j, reg_id, weight=cellss["reg_w_innode_list"][w_index])
        for j in range(cellss["regnr"]):
            reg_id2 = cellss["innr"] + j
            w_index = i * cellss["regnr"] + j
            netgraph.add_edge(reg_id2, reg_id, weight=cellss["reg_w_regnode_list"][w_index])

    for i in range(cellss["outnr"]):
        out_id = cellss["innr"] + cellss["regnr"] + i
        netgraph.add_node(out_id, type="out", subset=2, threshold=cellss["out_threshold_list"][i])
        for j in range(cellss["regnr"]):
            reg_id = cellss["innr"] + j
            w_index = i * cellss["regnr"] + j
            netgraph.add_edge(reg_id, out_id, weight=cellss["out_w_regnode_list"][w_index])

    return netgraph


if __name__ == "__main__":
    main()
