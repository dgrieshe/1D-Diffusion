#last edited by Miriam Rathbun on 6/18/2016
#This script plots all of the ouput from vector_practice separate runs on one figure
print('plotting all figures together >> figure_all.png')

#imports 
import sys, os
import matplotlib.pyplot as plt
import numpy as np

#variables from bash
n = int(sys.argv[1])
length = int(sys.argv[2])

with open("./output/output.text") as f:
	data = f.read()
	
data = data.split('\n')
plt.plot(data[0:n])
for i in range (1, (len(data)-1)/n):
	plt.plot(data[i*n:(i+1)*n])
plt.savefig('./output/figure_all')

print('done')


