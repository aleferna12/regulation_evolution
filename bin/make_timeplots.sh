runsdir=$1
outputdir=$2

for rundir in "$runsdir"/*
do
	python run_script.py plot_deaths "$rundir/cellgravedata_" "$outputdir/$(basename "$rundir")_death.html" -b 1000000 &
	python run_script.py plot_timeline "$rundir/celldata_" -o "$outputdir/$(basename "$rundir").html" -a "$outputdir/$(basename "$rundir")_adh.html" -n 1000 &
done
	