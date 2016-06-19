#!/bin/bash
#not actually sure what this first line does, but it is necessary
#last edited by Miriam Rathbun on 6/18/2016
#script reads an input file of known structure. It loops through all .inp files present 

for f in `ls *.inp`
do
   name=$f							#finds input files in the folder
   echo $name
 
   #finds the line where the word "length" appears, and prints column 2 of that line
   length=`awk '/length/{print $2}' $name`
   rightBC=`awk '/rightBC/{print $2}' $name` 
   leftBC=`awk '/leftBC/{print $2}' $name`
   numgroups=`awk '/numgroups/{print $2}' $name`
   numbins=`awk '/numbins/{print $2}' $name`
 
   i=1
   let j=i+1
   while [ $i -le $numbins ]; do
	   let j=i+1
	   let p=i-1
       xs=`awk '/xs/{print $"'"$j"'"}' $name`
       vector[$p]=$xs
       let i=i+1
   done

python ./diffusion.py $length $rightBC $leftBC $numgroups $f $numbins ${vector[*]}
done


#combines plots
python ./plotter.py $numbins $length