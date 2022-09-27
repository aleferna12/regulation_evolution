import sys
import os
import shutil
from pathlib import Path

EXEC_TARGET = "cell_evolution"
ROOT_PATH = Path(__file__).parent.parent.resolve()
SRC_PATH = ROOT_PATH / "src"
DATA_PATH = ROOT_PATH / "data"

assert len(sys.argv) == 4
rundir, nrruns = Path(ROOT_PATH / sys.argv[1]), sys.argv[2]
parfile = rundir / sys.argv[3]
execfile = rundir / EXEC_TARGET
assert not rundir.exists()

rundir.mkdir()
for file in list(SRC_PATH.iterdir()) + list(DATA_PATH.iterdir()):
    shutil.copy(file, rundir)

os.chdir(rundir)
os.system("qmake-qt4")
os.system("make")
os.system(f'"{execfile.resolve()}" "{parfile.resolve()}"')
