import json
import sys
import os
import logging
import re
from pathlib import Path
from ete3 import Tree, TreeStyle
from random import shuffle, seed
from parse_neighbours import read_neighbours
from colorir import *

_cs = Palette.load()


def main():
    logging.basicConfig(level=logging.INFO)
    season = int(sys.argv[1])
    treepath = Path(sys.argv[2]).resolve()
    neighpath = Path(sys.argv[3]).resolve()
    outpath = Path(sys.argv[4]).resolve()
    min_cluster = int(sys.argv[5]) if len(sys.argv) > 5 else 2

    logging.info("Reading tree file")
    tree = Tree(str(treepath))
    root_point = tree.get_midpoint_outgroup()
    tree.set_outgroup(root_point)

    logging.info("Reading neighbours file")
    if os.path.isfile(neighpath):
        with open(neighpath, "r") as file:
            neighs = json.load(file)[str(season)]
    else:
        logging.warning("Neighbours file not found")
        logging.info(f"Try to build neighbours from '{neighpath}'")
        neighs = read_neighbours(neighpath, season_filter=[season])[season]

    logging.info(f"Writing tree to '{outpath}'")
    plot_tree(tree, neighs, str(outpath), min_cluster)
    logging.info("Finished")


def get_cluster_colors(clusters, min_cluster):
    clusters = [c for c in clusters if len(c) >= min_cluster]
    palpath = Path(__file__).resolve().parent.parent.parent / "data"
    pal = StackPalette.load("carnival", palettes_dir=palpath)
    # Only allocate colors for clusters with at least min_cluster neighbours
    colors = PolarGrad(pal, color_sys=HCLuv).n_colors(len(clusters))
    # Shuffle colors to avoid confusion of similar colors
    seed(0)
    shuffle(colors)
    return colors


def plot_tree(tree: Tree, clusters, outpath, min_cluster=2):
    colors = get_cluster_colors(clusters, min_cluster)
    leaf_color = {}
    for cluster, color in zip(clusters, colors):
        for leaf in cluster:
            leaf_color[str(leaf)] = color

    for node in tree.traverse():
        node.img_style["size"] = 0
        if node.is_leaf():
            # rapidnj adds these for no reason
            node.name = node.name.strip("\'")
            node.img_style["bgcolor"] = leaf_color.get(node.name, _cs.white)

    ts = TreeStyle()
    ts.mode = 'c'
    tree.render(outpath, tree_style=ts, w=1000)


def figtree_nexus_str(newick, clusters, min_cluster=2):
    """Creates a color-coded nexus string that can be parsed by FigTree.

    This approach was abandoned because the clusters couldn't be distinguished in the end result.
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
