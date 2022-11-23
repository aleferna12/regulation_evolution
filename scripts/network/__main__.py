import argparse
import re
import logging
import sys
from functools import partial
from pathlib import Path
from parse_ancestry import parse_ancestry
from plot_tree import plot_tree


def main():
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Create tree plots from a directory with runs "
                                                 "(that must follow a directory naming "
                                                 "convention).")
    parser.add_argument("--runsdir")
    parser.add_argument("--outdir")
    args = parser.parse_args()

    runsdir = Path(args.runsdir)
    outdir = Path(args.outdir)
    if not runsdir.is_dir():
        logging.warning("'runsdir' is not valid directory")
        logging.warning("Aborting")
    if not outdir.is_dir():
        logging.warning("'outdir' is not valid directory")
        logging.warning("Aborting")

    for rundir in runsdir.iterdir():
        logging.info(f"Creating trees for run: '{rundir.name}'")
        netdir = rundir / "networkdir_"
        latest_time = 0
        for anc_file in netdir.glob("anc*"):
            s_time = int(re.search(r"(?<=t)\d+", anc_file.name).group())
            latest_time = max(latest_time, s_time)

        outrundir = outdir / rundir.name
        outrundir.mkdir()
        (outrundir / "trees").mkdir()
        (outrundir / "trees_last").mkdir()

        trees = parse_ancestry(netdir, outrundir / "trees", True, True, True, True, 5)
        trees_last = parse_ancestry(netdir, outrundir / "trees_last", False, True, True, True, 5)
        plot_foo = partial(
            plot_tree,
            timepoint=latest_time,
            neighpath=netdir,
            min_cluster=2,
            reroot=False,
            colored=True
        )
        for tree in trees:
            plot_foo(
                treepath=outrundir / "trees" / (tree + ".newick"),
                outpath=outrundir / (tree + ".svg"),
            )
        for tree in trees_last:
            plot_foo(
                treepath=outrundir / "trees_last" / (tree + ".newick"),
                outpath=outrundir / f"{tree}_last.svg"
            )


if __name__ == "__main__":
    main()