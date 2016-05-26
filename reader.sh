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


python ./vector_practice.py $length $rightBC $leftBC $numgroups $f
done

#combines plots
python ./plotter.py $length