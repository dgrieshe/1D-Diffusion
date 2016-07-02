#last edited by Miriam Rathbun on 7/2/2016
#this script is run from reader.sh. It solves the discrete diffusion equations
#this script can run multiple input files with restrictions



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from nuclide import *
from diffOpts import *



#variables
#name of the input file
f = sys.argv[1]		
options=DiffusionOpts1D()
options.read(f)
#retains only the number associated with the input file. Ex: "input1" becomes "1" 
#This was made to make the output name nice
name = f[len(f)-5]										 
#cross sections for each bin.
xs=np.zeros(options.numBins)
#diffusion coefficients. Should be read from input file, eventually. D=1/(3*xs)
diff_coef=np.zeros(options.numBins+1)
#Below, you could specify a different nuclide for given values of i, if desired
#for example, if i<options.length/2: M=Nuclide('U235') else: M=Nuclide('U238')
for i in range(0,options.numBins):
	M= Nuclide('U235')
	M.read()
	xs[i]= M.absXs
	diff_coef[i]=1/(3*M.absXs)
#arbitrary source term, constant            
source=[5]
                                                                                  

###############################################################################
#Filling matrices A (linear system of diffusion equations) and B (source)
A=np.zeros((options.numBins,options.numBins)) 										
B=np.zeros((options.numBins,1))										
det=np.linalg.det(A) 
#if det=0, A is singular, signal an error, eventually.
while det==0:											
	
   for row in range(0,options.numBins):
	   for col in range(0,options.numBins):
	   	   if row == col:
	   	   	   if row == 0:
	   	   	   	   #this uses the ghost bin. I don't know what to assign it yet, so I'm giving it a value of 1. 
	   	   	   	   A[row,col]=2*diff_coef[row]*1/(options.delta*1/+options.delta*diff_coef[row])+xs[row]*options.delta
	   	   	   else:
	   	   	   	   if row == options.numBins-1:
	   	   	   	   	   #this uses the ghost bin. I don't know what to assign it yet, so I'm giving it a value of 1.
	   	   	   	   	   A[row,col]=2*diff_coef[row]*1/(options.delta*1/+options.delta*diff_coef[row])+xs[row]*options.delta
	   	   	   	   else:
	   	   	   	       A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])+xs[row]*options.delta
	   	   else:                                                            
	   	   	   if col == row-1:
	   	   	   	   A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
	   	   	   else:
	   	   	   	   if col == row+1:
	   	   	   	   	   A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(options.delta*diff_coef[col]+options.delta*diff_coef[col-1]) 
	   B[row,0]=source[0]*options.delta		    
   det=np.linalg.det(A)
###############################################################################   
 
 
###############################################################################
#calculating the solution x to Ax=B
Ainv=np.linalg.inv(A)
x=np.dot(Ainv,B)
###############################################################################



###############################################################################
#plotting solution x
plt.plot(x)
plt.ylabel('some numbers')
plt.savefig('./output/figure'+name)

if name=='1':
    f = open('./output/output.text', 'w')
    f.close()

f = open('./output/output.text', 'a')

i=0
while i < options.numBins:
   print >> f, float(x[i])
   i = i+1
f.close()
###############################################################################

                   