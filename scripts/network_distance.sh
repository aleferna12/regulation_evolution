#!/bin/bash

orifile=$1

ser=$(ls  | tail -10 |head -1)
timepoint=${ser%_*}

scriptdir="/home/renske/origins_postdoc/multicell/regulation_evolution/scripts"

outfile=distances_$timepoint.txt
echo -n "" >$outfile

#declare -a filemontage
#loop through concentrations and divisions
#divs=0
sumsum=0
stdev=0.
numnum=0
 
#echo "processing $orifile"
count=0
for fi in ${timepoint}_c*.txt; do
 if [ "$count" -gt "0" ]; then
   diff=$(/usr/bin/python3 ${scriptdir}/count_2_networks_dist.py $orifile $fi)
   echo $diff >> $outfile
   sumsum=$(bc -l <<<"$sumsum+$diff")
   stdev=$(bc -l <<<"$stdev+$diff*$diff")
   ((numnum++))
 fi
 ((count++))
done

sumsum=$(bc -l <<<"$sumsum/$numnum") #the average
stdev=$(bc -l <<<"$stdev/$numnum-$sumsum*$sumsum") #variance
echo "$sumsum $stdev"
