#last edited by Miriam Rathbun on 7/18/2016
#This script solves the discrete diffusion equations
#this script can run multiple input files with restrictions



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np    
from nuclide import *
from diffOpts import *
from plotter import * 
from GetFileName import *

fn=FileName()
fn.GetFileName()

for f in fn.listfn:
	print f
	#retains only the number associated with the input file. Ex: "input1" becomes "1" 
	#This was made to make the output name nice
	name = f[len(f)-5]
	
	#variables
	u_g=0.33
	
	options=DiffusionOpts1D()
	options.read(f)
	
	#cross sections for each bin.
	xs=[]
	diff_coef=[]
	source=np.zeros(options.numBins*options.numGroups)
	scat=np.zeros((options.numGroups,options.numGroups))
	
	#creates slab that's all U-235
	for k in range(1,options.numGroups+1):
		M= ('M%i' %k)
		M=Nuclide('U235(%i)' %k)
		M.read()
		for i in range(0,options.numBins):
			xs.append(M.absXs)
			diff_coef.append(1/(3*(M.totalXs-u_g*M.scatXs)))
		for j in range(0,options.numGroups):
			scat[j,k-1]=eval('M.scatXs%i' %(j+1))
			
			
				
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
		
	source[10]=0.1
	#for i in range(0,len(source)):
	#	source[i]=0.1
	
###############################################################################
#Filling matrices A (linear system of diffusion equations) and B (source)

	
	A=np.zeros((options.numBins*options.numGroups,options.numBins*options.numGroups))
	B=np.zeros((options.numBins*options.numGroups,1))
	

	for k in range(1,options.numGroups+1):
		for row in range(options.numBins*(k-1),options.numBins*k):
			for col in range(options.numBins*(k-1),options.numBins*k):
				if row == col:
					if row == options.numBins*(k-1):
						A[row,col]=2*diff_coef[row]/(options.delta*1+4*diff_coef[row])+xs[row]*options.delta #+2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])
					elif row == options.numBins*k-1:
						A[row,col]=2*diff_coef[row]/(options.delta*1+4*diff_coef[row])+xs[row]*options.delta #+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
					else:
						A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])+xs[row]*options.delta
				elif col == row-1:
					A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
				elif col == row+1:
					A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(options.delta*diff_coef[col]+options.delta*diff_coef[col-1]) 
		for row in range(options.numBins,options.numBins*options.numGroups):
			for col in range(0,options.numBins*options.numGroups):
				if row == col+options.numBins*k:
					A[row,col]=-scat[k,k-1]
	for row in range(0,len(source)):
		B[row,0]=source[row]*options.delta

###############################################################################   

    
     
###############################################################################
#calculating the solution x to Ax=B   
	Ainv=np.linalg.inv(A)
	x=np.dot(Ainv,B)
	print x[0]

	
	results=Plotter()
	results.plot(x,1,options.numBins,options.numGroups,name)
    
###############################################################################
    
    
               