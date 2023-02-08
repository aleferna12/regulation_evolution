#!/bin/bash

dire=$1
start=$2		#typical value when finished: 9 000 0000
dire=${dire%/*}
cd $dire

echo $dire

mkdir formovie/


#copy last few seasons to new dir
for i in tau*.png; do 

  nr=$(echo "$i" | tr -dc '0-9') 
  nr2=${nr##+(0)} 
  if [ "$nr2" -ge "$start" ]; then 
    #echo "$nr2" 
    cp $i formovie/ 
  fi 

done

cd formovie/
#make actual movie

ffmpeg -framerate 16 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p $dire.mp4

mv $dire.mp4 ../

cd ..

rm -rf formovie/

