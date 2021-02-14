#!/bin/bash

#very basic script to run x simulations with a particular parameter setting


#echo "Make sure the parfile is ok for running 1 season, no pictures"

networkfile=$1
parfile=$2
nrcores=$3

exe="./cell_evolution"
exewrap="basicwrap.sh"

output_sigmas () {
  while IFS=$'\t' read -r -a line
  do

    echo ${line[0]}
    
  done < $1
}

#turn this into a gnu parallel?
parallel --gnu -j $nrcores $exewrap $exe $parfile <<< "$(seq $networkfile)"


