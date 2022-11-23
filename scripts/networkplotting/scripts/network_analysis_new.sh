#!/bin/bash

netfile=$1


## parameters to pass as input to the genomecode for running
## note: I copied the executable of the genome code to my scriptdir, out of laziness
MINGRAD=0
MAXGRAD=107  # roof(par.gradscale * sqrt(2) * latt_side / 100)
GRADSTEP=4
MINFOOD=0
MAXFOODIVS=100
FOODDIVSSTEP=10

grad=$MINGRAD
food=$MINFOOD
scriptdir=/home/aleferna/CProjects/Projects/regulation_evolution/scripts/networkplotting/scripts

echo "processing $netfile"
#convert -size 20x20 canvas:firebrick red.png
#convert -size 20x20 canvas:royalblue blue.png
#: <<'COMMENT'
#store the states of the network

outfile=states_${netfile%.*}.dat
echo -n "" >"$outfile"

#declare -a filemontage
#loop through concentrations and food

printf "0">>"$outfile"
while [ "$food" -le "$MAXFOODIVS" ];
do
  printf "\t%d" $food >> "$outfile"
  food=$(bc <<< "$food+$FOODDIVSSTEP")
done
printf "\n" >> "$outfile"
food=$MINFOOD

while [ "$grad" -le "$MAXGRAD" ]; do 
  
  printf "%d" $grad >>"$outfile"
  while [ "$food" -le "$MAXFOODIVS" ]; do

   #run_network returns the network output over 100 time steps
   mapfile -t my_array < <( ${scriptdir}/run_network 572843 "$netfile" $grad $food )
   sumsum=0
   for el in "${my_array[@]:80:20}"; do
     #echo $el
     sumsum=$(bc<<<"$sumsum+${el:0:1}")
   done
   #echo $sumsum
   if [ "$sumsum" -lt "19" ]; then
     state=0
     #echo "It is migrating, $sumsum"
     #filemontage+=("red.png")
   else
     state=1
     #filemontage+=("blue.png")
     #((divs++))
   fi
   printf "\t%d" $state >>"$outfile"
   food=$(bc <<< "$food+$FOODDIVSSTEP")
  done
  printf "\n">>"$outfile"

grad=$(bc <<< "$grad+$GRADSTEP")
food=$MINFOOD
done

#xval=$(bc <<< "($MAXDIV+1)/$DIVSTEP")
#echo ${filemontage[@]}
#montage -geometry +0+0 -background '#000000' -tile ${xval}x ${filemontage[@]} tempimage.png
#divnr=$(printf "%04d\n" $divs)
#convert -flip tempimage.png div_${divnr}_${netfile%.*}_im.png 
#rm -rf tempimage.png
#COMMENT

#python ${scriptdir}/plot_netanalysis.py $outfile ${netfile%.*} $MINGRAD
### call this after all individuals have been done
#python ${scriptdir}/plot_netanalysis_jan.py $outfile ${netfile%.*} 1 $MINGRAD

