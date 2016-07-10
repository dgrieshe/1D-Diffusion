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
#j is the group number
j=options.numGroups
#retains only the number associated with the input file. Ex: "input1" becomes "1" 
#This was made to make the output name nice
name = f[len(f)-5]		
                     
#cross sections for each bin.
xs=np.zeros(options.numBins*j)
diff_coef=np.zeros(options.numBins*j)      
                                          

    
#creates slab half U-235 and half U-238
for i in range(0,options.numBins*j):
	M1=Nuclide('U235(1)') 
	M2=Nuclide('U235(2)')                 
	M1.read()
	M2.read()                                
	if i<options.numBins:
		xs[i]=M1.absXs
		diff_coef[i]=1/(3*M1.absXs)
	else:
		xs[i]=M2.absXs
		diff_coef[i]=1/(3*M2.absXs)

		
			
			
	
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
	
for row in range(0,options.numBins*j):
 for col in range(0,options.numBins*j):
 	   if row == col:
 	   	   if row == 0 or row == options.numBins or row == options.numBins-1 or row == options.numBins*j-1:
 	   	   	   #this uses the ghost bin with delta=4 and diff_coef=1 
 	   	   	   A[row,col]=2*diff_coef[row]*1/(options.delta*1/+4*diff_coef[row])+xs[row]*options.delta
 	   	   else:
 	   	   	   A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])+xs[row]*options.delta
 	   else:                                                            
 	   	   if col == row-1 and row != options.numBins:
 	   	   	   A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
 	   	   else:
 	   	   	   if col == row+1 and col != options.numBins:
 	   	   	   	   A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(options.delta*diff_coef[col]+options.delta*diff_coef[col-1]) 
 B[row,0]=source[0]*options.delta	
###############################################################################   


 
###############################################################################
#calculating the solution x to Ax=B
Ainv=np.linalg.inv(A)
x=np.dot(Ainv,B)
###############################################################################
  


###############################################################################
#plotting solution x
plt.plot(x)
plt.ylabel('flux')
plt.savefig('./output/figure'+name)

if name=='1' and j==1:
    f = open('./output/output.text', 'w')
    f.close()

f = open('./output/output.text', 'a')

i=0
while i < options.numBins:
   print >> f, float(x[i])
   i = i+1
f.close()
###############################################################################
               