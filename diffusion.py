#last edited by Miriam Rathbun on 7/9/2016
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
#j is the number of groups
j=options.numGroups
#retains only the number associated with the input file. Ex: "input1" becomes "1" 
#This was made to make the output name nice
name = f[len(f)-5]		


#cross sections for each bin.
xs= []
diff_coef=[]     
                                          

    
#creates slab that's all U-235
for k in range(1,j+1):
	M= ('M%i' %k)
	M=Nuclide('U235(%i)' %k)
	M.read()
	for i in range(0,options.numBins):
		xs.append(M.absXs)
		diff_coef.append(1/(3*M.absXs))
			
		
	
	#if i<=options.numBins/2:
	#	M=Nuclide('U235(%i)' % j)
	#	M.read()
	#	xs[i]= M.absXs                                  
	#	diff_coef[i]=1/(3*M.absXs)
	#else:
	#	M= Nuclide('U238(%i)' %j)
	#	M.read()
	#	xs[i]= M.absXs                                       
	#	diff_coef[i]=1/(3*M.absXs)       
		
#arbitrary source constant source term           
source=[5]                                                                          

###############################################################################
#Filling matrices A (linear system of diffusion equations) and B (source)
A=np.zeros((options.numBins*j,options.numBins*j)) 										
B=np.zeros((options.numBins*j,1))																				
	
for k in range(1,j+1):
		for row in range(options.numBins*(k-1),options.numBins*k):
			for col in range(options.numBins*(k-1),options.numBins*k):
				if row == col:
					if row == options.numBins*(k-1) or row == options.numBins*k-1:
						A[row,col]=2*diff_coef[row]*1/(options.delta*1/+4*diff_coef[row])+xs[row]*options.delta
					else:
						A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])+xs[row]*options.delta
				else:
					if col == row-1:
						A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
					else:
						if col == row+1:
							A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(options.delta*diff_coef[col]+options.delta*diff_coef[col-1]) 
			B[row,0]=source[0]*options.delta	
###############################################################################   


 
###############################################################################
#calculating the solution x to Ax=B
Ainv=np.linalg.inv(A)
x=np.dot(Ainv,B)

for k in range(1,j+1):
	l=0
	x_group= np.zeros(options.numBins)
	for i in range(options.numBins*(k-1),options.numBins*k):
		x_group[l]=x[i]
		l=l+1
	plt.plot(x_group)
	plt.ylabel('flux')
	plt.savefig('./output/figure'+name)

###############################################################################


               