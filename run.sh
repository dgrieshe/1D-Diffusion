#!/bin/bash
#last edited by Miriam Rathbun on 8/3/2016
#runs diffusion.py n times

n=100
i=0
while [ $i -lt $n ]
do
	python ./diffusion.py
	let i=i+1
done