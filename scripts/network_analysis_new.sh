#!/bin/bash

netfile=$1

MINGRAD=0
MAXGRAD=28
GRADSTEP=1
MINDIV=0
MAXDIV=3
DIVSTEP=1

grad=$MINGRAD
div=$MINDIV
scriptdir=~/origins_postdoc/multicell/regulation_evolution/scripts/

echo "processing $netfile"
#convert -size 20x20 canvas:firebrick red.png
#convert -size 20x20 canvas:royalblue blue.png
#: <<'COMMENT'
#store the states of the network

outfile=states_${netfile%.*}.dat
echo -n "" >$outfile

#declare -a filemontage
#loop through concentrations and divisions
#divs=0

printf "0\t0\t1\t2\t3\n">>$outfile
while [ "$grad" -le "$MAXGRAD" ]; do 
  
  printf "%d" $grad >>$outfile
  while [ "$div" -le "$MAXDIV" ]; do 

   #run_network returns the network output over 100 time steps
   mapfile -t my_array < <( ${scriptdir}/run_network 572843 $netfile $grad $div )
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
   printf "\t%d" $state >>$outfile
   div=$(bc <<< "$div+$DIVSTEP")
  done
  printf "\n">>$outfile

grad=$(bc <<< "$grad+$GRADSTEP")
div=$MINDIV
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

