import sys
import logging
import re
from pathlib import Path
from ete3 import Tree, TreeStyle, AttrFace
from random import shuffle, seed
from fileio import parse_cell_data
from plot_timeline import get_adhering_clusters, get_cluster_colors, CellCluster
from colorir import *

logger = logging.getLogger(__name__)
_cs = Palette.load()


def main():
    logging.basicConfig(level=logging.INFO)

    treepath = sys.argv[1]
    datafile = sys.argv[2]
    outpath = sys.argv[3]
    min_cluster = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    colored = False if len(sys.argv) > 5 and sys.argv[5] in ["false", '0'] else True
    reroot = True if len(sys.argv) > 6 and sys.argv[6] in ["true", '1'] else False

    logger.info("Reading tree file")
    tree = Tree(treepath, format=5)
    # When tree is constructed by neighbour joining they are unrooted so we must estimate
    if reroot:
        mid = tree.get_midpoint_outgroup()
        tree.set_outgroup(mid)

    celldf = parse_cell_data(datafile)
    clusters = get_adhering_clusters(celldf)
    plot_tree(tree, clusters, outpath, min_cluster, colored)

    logger.info("Finished")


def plot_tree(tree: Tree, clusters: CellCluster, outpath, min_cluster=2, colored=True):
    logger.info(f"Plotting tree to '{outpath}'")

    colors = get_cluster_colors(clusters, min_cluster)
    leaf_color = {}
    for cluster, color in zip(clusters, colors):
        for leaf in cluster:
            leaf_color[str(leaf)] = color

    last_timepoint = getattr(tree.get_farthest_leaf()[0], "time", None)
    for node in tree.traverse():
        node.img_style["size"] = 0
        if node.is_leaf():
            # rapidnj adds these for no reason
            node.name = node.name.strip("\'")
            face = AttrFace("name")
            node.add_face(face, column=0)
            # Color unicellular nodes dark gray and doesn't color dead-end nodes
            if colored and getattr(node, "time", None) == last_timepoint:
                if node.name not in leaf_color:
                    raise ValueError("'clusters' doesn't contain one or more terminal leaves "
                                     "from 'tree', check if clusters and tree come from the same "
                                     f"time point: {last_timepoint}")
                node.img_style["bgcolor"] = leaf_color[node.name]
                node.img_style["hz_line_color"] = leaf_color[node.name]

    ts = TreeStyle()
    ts.root_opening_factor = 0.1
    ts.show_leaf_name = False
    ts.mode = 'c'
    tree.render(str(outpath), tree_style=ts, w=250, units="mm")


def figtree_nexus_str(newick, clusters, min_cluster=2):
    """Creates a color-coded nexus string that can be parsed by FigTree.

    This approach was abandoned because the clusters couldn't be distinguished in the end result.
    That said, the distances between the branches are much easier to see in FigTree thanks to the
    "radial" view (but you don't need this function for that, just open the normal tree in FigTree
    or Dendroscope).
    """
    def sub_color(match):
        number = match.group(1)
        if number not in leaf_color:
            return match.group()
        lcolor = leaf_color.get(number, _cs.black)
        return f"{number}[&!color={lcolor}]:"

    colors = get_cluster_colors(clusters, min_cluster)
    leaf_color = {}
    for cluster, color in zip(clusters, colors):
        for leaf in cluster:
            leaf_color[str(leaf)] = color
    figtree_fmt = re.sub(r"(\d+):", sub_color, newick)

    return (
        "#NEXUS\n"
        "Begin Trees;\n"
        f" Tree tree1 = [&R] {figtree_fmt}\n"
        "End;\n"
    )


if __name__ == "__main__":
    main()