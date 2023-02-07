#!/bin/bash

# Run single individual
#
# Usage: bash exewrap.sh EXE PARFILE PARENTDIR CHILDDIR

# module purge
# module load ...

set -e



exe=$1                          #executable; elli in my case
parfile=$2                      #parfile with mutation rates 
rundir=$3			#directory in which to put stuff
runnr="$4"                        #the copy to be made


$exe $parfile -datafile $rundir/data_$runnr.txt -datadir $rundir/run_$runnr -backupdir $rundir/bckp_$runnr 
#echo "run is done"

