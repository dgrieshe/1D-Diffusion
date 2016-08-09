#last edited by Miriam Rathbun on 8/3/2016
#This script solves the discrete diffusion equations



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np  
from nuclide import *
from diffOpts import *
from plotter import * 
from GetFileName import *
from solver import *
from Construct import *

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
	sourceGroup=[]
	
###############################################################################	
	
	#creates slab that's all U-235
	#loops through the groups
	#this section also fills other items which depend on groups, see below
	for k in range(1,options.numGroups+1):
		M= ('M%i' %k)
		M=Nuclide('U235(%i)' %k)
		M.read()                    
		for i in range(0,options.numBins):
			xs.append(M.totalXs)
			diff_coef.append(1/(3*(M.totalXs-u_g*M.scatXs)))
			
		#fills the transition matrix/scattering kernel
		for j in range(0,options.numGroups):  
			scat[j,k-1]=eval('M.scatXs%i' %(j+1))
			
			
###############################################################################

		#previous code to make a half-slab of U235 and U238
				
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
    	
	
###############################################################################
#Builds the linear system of equations
			
	A=Construct()
	A.constructA(options, diff_coef, scat, xs)
	
#sets the initial value of the source
	for i in range (0,len(source)):
		source[i]=1/options.length
		
		
###############################################################################
#calculates the solution x to Ax=B 	
	
	
	A.invertA()
	sol=Solve()
	sol.solve(options, A.inv, source)
	
	results=Plotter()
	results.plot(sol.x,1,options.numBins,options.numGroups,name)
    
###############################################################################
