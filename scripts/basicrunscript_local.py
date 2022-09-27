import sys
import os
import shutil
from datetime import datetime
from pathlib import Path

EXEC_TARGET = "cell_evolution"
ROOT_PATH = Path(__file__).parent.parent.resolve()
SRC_PATH = ROOT_PATH / "src"
DATA_PATH = ROOT_PATH / "data"

assert len(sys.argv) == 4
rundir, nrruns = Path(ROOT_PATH / sys.argv[1]), sys.argv[2]
parfile = rundir / sys.argv[3]
execfile = rundir / EXEC_TARGET

if not rundir.exists():
    rundir.mkdir()
for file in list(SRC_PATH.iterdir()) + list(DATA_PATH.iterdir()):
    old = rundir / file.name
    if old.exists():
        os.remove(old)
    shutil.copy(file, rundir)

os.chdir(rundir)
os.system("qmake-qt4")
os.system("make")
os.system(f'"{execfile.resolve()}" "{parfile.resolve()}"')

now = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
for path in ["backup", "data_film2", "data_cellcount.txt"]:
    os.rename(rundir / path, f"{str(rundir / path)}_{now}")
