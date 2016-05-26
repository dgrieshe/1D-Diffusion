#last edited by Miriam Rathbun on 5/26/2016
#This script plots all of the ouput from vector_practice separate runs on one figure
print('plotting all figures together >> figure_all.png')

#imports 
import sys, os
import matplotlib.pyplot as plt
import numpy as np

#variables from bash
length = int(sys.argv[1])

with open("./output/output.text") as f:
	data = f.read()
	
data = data.split('\n')
plt.plot(data[0:length])
for i in range (1, (len(data)-1)/length):
	plt.plot(data[i*length:(i+1)*length])
plt.savefig('./output/figure_all')

print('done')


