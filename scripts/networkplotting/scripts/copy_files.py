import random
import re
import shutil
import sys
from pathlib import Path


def main():
    sourcedir = Path(sys.argv[1]).resolve()
    outputdir = Path(sys.argv[2]).resolve()
    timepoint = int(sys.argv[3])
    k = int(sys.argv[4])

    copy_files(sourcedir, outputdir, timepoint, k)


def copy_files(sourcedir, targetdir, timepoint, k):
    if int(timepoint) == -1:
        timepoints = set(int(re.search(r"\d+", filepath.name).group())
                         for filepath in Path(sourcedir).iterdir())
        timepoint = sorted(timepoints)[-1]

    files = [filepath for filepath in Path(sourcedir).iterdir()
             if f"{timepoint:0>9}" in filepath.name and "_c" in filepath.name]
    for file in random.choices(files, k=k):
        shutil.copy(file, targetdir)


if __name__ == "__main__":
    main()
