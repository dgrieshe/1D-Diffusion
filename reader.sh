#!/bin/bash
#not actually sure what this first line does, but it is necessary
#last edited by Miriam Rathbun on 7/2/2016
#Loops through all .inp files present in the folder

for f in `ls *.inp`
#finds input files in the folder
do
   name=$f
   echo $name
python ./diffusion.py $f
done


#combines plots
python ./plotter.py $f