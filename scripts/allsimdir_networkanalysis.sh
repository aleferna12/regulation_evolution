#!/bin/bash

#this script assesses the phenotype for all simulations in the dir
season=$1
redo=$2

scriptdir=~/origins_postdoc/multicell/regulation_evolution/scripts/
for dir in movie_sim*; do 
    cd ${dir}/networks/

    echo $dir
    #find last timepoint
    ser=$(ls  | tail -10 |head -1)
    timepoint=${ser%_*}
    
    done=$(ls states_t$timepoint* | wc -l)
    tau=$(ls t$timepoint* | wc -l)
    echo "assessed $done out of $tau"
    
    if [ "$redo" -gt "0" || "$done" -lt "$tau" ]; then
        rm -rf states_t$timepoint*
    else        
        cd -
        continue
    fi

    #go through all individuals in the last pop, make matrix
    for net in t${timepoint}_c*.txt; do 
        echo $net 
        if [ "$count" -gt "0" ]; then 
            bash ~/origins_postdoc/multicell/regulation_evolution/scripts/network_analysis_new.sh ${net} 
        fi
        ((count++))
    done
    
    #print all matrices, find consensus
    simname=${dir#movie_*}
    simname=${simname%/*}
    python3 ~/origins_postdoc/multicell/regulation_evolution/scripts//plot_netanalysis_jan.py ${simname} 1 states_t$timepoint*.dat > consensusdata.txt
    
    #gather all images
    all_ims=(); for im in div_*pdf; do all_ims+=( "$im" ); ((count++)); done
    montage -geometry +4+4 -background '#000000' -tile 50x ${all_ims[@]} allmatrixplots.pdf

    cp consensus* ../../

    cd -
done



