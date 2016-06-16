#!/bin/bash
#not actually sure what this first line does
#last edited by Miriam Rathbun on 5/26/2016
#script reads an input file of known structure. It loops through all .inp files in the folder and runs the python script for each one of those. 

for f in `ls *.inp`
do
   name=$f							#finds input files in the folder
   echo $name

   length=`awk '/length/{print $2}' $name` 					#finds the line where the word "length" appears, and prints column 2 of that line
   rightBC=`awk '/rightBC/{print $2}' $name` 
   leftBC=`awk '/leftBC/{print $2}' $name`
   numgroups=`awk '/numgroups/{print $2}' $name`
   numbins=`awk '/numbins/{print $2}' $name`
   xs1=`awk '/xs/{print $2}' $name`
   xs2=`awk '/xs/{print $3}' $name`
   xs3=`awk '/xs/{print $4}' $name`
   xs4=`awk '/xs/{print $5}' $name`


python ./diffusion.py $length $rightBC $leftBC $numgroups $f $numbins $xs1 $xs2 $xs3 $xs4 $xs5
done

#combines plots
python ./plotter.py $numbins