#!/bin/bash

#this script is derived from allsimdir_networkanalysis.sh
#it makes the network plot for all individuals at a particular timepoint; if you specify timepoint as -1, then it will pick the last one
# it calls another script to actually run the genome dynamics: network_analysis_new.sh
#and then it calls a pythonscript to make all the pictures and extract a consensus network: plot_netanalysis_jan.py
dir=$1
redo=$2
timepoint=$3
samples=$4
copy_from=$5
# usage: ./onesimdir... ../data 0 -1 100 <.../networkdir_>

scriptdir=/home/aleferna/CProjects/Projects/regulation_evolution/scripts/networkplotting/scripts
#for dir in movie_sim*; do
    cd ${dir}/
    if [ $# -gt 3 ]
      then
        python3 "$scriptdir/copy_files.py" "$copy_from" . "$timepoint" "$samples"
      fi

    echo "$dir"
    #find last timepoint
    if [ "${timepoint}" -le "0" ]; then
      ser=$(ls  | tail -10 |head -1)
      timepoint=${ser%_*}
      echo $timepoint
    else
      timepoint=$(printf "t%010d\n" ${timepoint})
    fi
    
    done=$(ls states_$timepoint* | wc -l)
    tau=$(ls $timepoint* | wc -l)
    echo "assessed $done out of $tau"
    
    if [ "$redo" -gt "0" ] || [ "$done" -lt "$tau" ]; then
        rm -rf states_$timepoint*
        echo "yeah"
    else        
        cd -
        continue
    fi
    count=0
    pwd
    #go through all individuals in the last pop, make matrix
    for net in ${timepoint}_c*.txt; do
        echo $net
        if [ "$count" -gt "0" ] && [ ! -f "states_${net%.*}.dat" ] ; then
            bash ${scriptdir}/network_analysis_new.sh ${net}
        fi
        ((count++))
    done
    
    #print all matrices, find consensus
    simname=${dir#movie_*}
    simname=${simname%/*}
    python3 ${scriptdir}/plot_netanalysis_jan.py ${simname} 1 states_$timepoint*.dat > consensusdata.txt
    
    #gather all images into one !NOT WORKING!
    all_ims=(); for im in div_*pdf; do all_ims+=( "$im" ); ((count++)); done
    montage -geometry +4+4 -background '#000000' -tile 50x "${all_ims[@]}" allmatrixplots.pdf

    # cp consensus* ../../

    cd -
#done



