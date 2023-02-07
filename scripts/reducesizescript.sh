#!/bin/bash

#this script will clean up your directories by removing most backup files and making a tarball of the datafiles and the pngs.
#these can be later be retrieved one by one or all at the same time
#Example: suppose you want file etc/apt/sources.list from etc.tar:
#tar -xf etc.tar etc/apt/sources.list
#Will extract sources.list and create directories etc/apt under the current directory.

dire=$1

dire=${dire%/*}
cd $dire

echo $dire

#remove superfluous backup files
for dir in backup_sim*; do

  cd $dir
  lastfile="$(ls | tail -1)"
  nr=$(echo "$lastfile" | tr -dc '0-9')
  nr2=$(bc <<< "$nr -5000000")

  for i in backup_*; do 
    NUMBER=$(echo "$i" | tr -dc '0-9') 

    if [ "$NUMBER" -lt "$nr2" ]; then 
      rm $i
    fi
  done
  cd ..
done

#tar all data files
echo "we are going to tar the datafiles"
tarred=$(ls data_sim*.txt | wc -l)
if [ "$tarred" -gt "0" ]; then
  #rm -rf datafiles.tgz
  tar -zcf datafiles.tgz data_sim*.txt
  wait
  rm -rf data_sim*.txt
fi

#make videos of last few seasons and tar the pngs and network files
for dir in movie_sim*; do

  echo "now in this dir: $dir"
  cd $dir

  video=$(ls *.mp4 | wc -l)

  #we still have to make a video
  if [ "$video" -eq "0" ]; then
    lastfile="$(ls | tail -1)"
    nr=$(echo "$lastfile" | tr -dc '0-9')
    nr2=$(bc <<< "$nr -5000000")
    nr3=$(bc <<< "$nr2/10000000*10000000") #filters the smaller numbers, leaving a nice round number. Assumes the simulation got to at least 10^7 steps.
    cd ..
    bash ~/origins_postdoc/multicell/regulation_evolution/scripts/moviescript.sh $dir $nr3
  cd $dir
  fi
  
  #now tar images
  tarred=$(ls tau*.png | wc -l)
  if [ "$tarred" -gt "0" ]; then
    #rm -rf datafiles.tgz
    tar -zcf images.tgz tau*.png
    wait
    rm -rf tau*.png
  fi
  
  
  #now tar networks
  if [ -d networks ]; then
    tar -zcf networks.tgz networks/
    wait
    rm -rf networks/
  fi
  cd ..
  
done

