#last edited by Miriam Rathbun on 8/15/2016
#This script solves the discrete diffusion equations



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np 
from diffOpts import *
from plotter import * 
from GetFileName import *
from solver import *
from Construct import *
from nuclide import *
from material import *

fn=FileName()
fn.GetFileName()

for f in fn.listfn:
	print f
	#retains only the number associated with the input file. Ex: "input1" becomes "1" 
	#This was made to make the output name nice
	name = f[len(f)-5]
	
	#read input and calculate macroscopic cross sections
	options=DiffusionOpts1D()
	options.read(f)
	
	M=Material()
	M.read()
	M.calc_macro()
	
	#variables
	u_g=0.33
	
	totXS=[]
	absXS=[]
	scatXS=[]
	fisXS=[]
	diffcoef=[]
	
	
	numDensity=np.zeros((options.numBins,3))
	numDensity[:,0]=M.NDfuel
	numDensity[:,1]=M.NDmod
	numDensity[:,2]=M.NDpoison
	source=np.zeros(options.numBins*options.numGroups)
	scat=np.zeros((options.numGroups,options.numGroups))
	sourceGroup=[]
	
	
###############################################################################	
	
#create slab and assign each bin its cross section by energy group

	for k in range(1,options.numGroups+1):
		for i in range(0,options.numBins):
			totXS.append(M.data['fuel']['totXS'][k]+M.data['moderator']['totXS'][k]+M.data['poison']['totXS'][k])
			scatXS.append(M.data['fuel']['scatXS'][k]+M.data['moderator']['scatXS'][k]+M.data['poison']['scatXS'][k])
			absXS.append(M.data['fuel']['absXS'][k]+M.data['moderator']['absXS'][k]+M.data['poison']['absXS'][k])
			fisXS.append(M.data['fuel']['fisXS'][k]+M.data['moderator']['fisXS'][k]+M.data['poison']['fisXS'][k])
			
			diffcoef.append(1/(3*(totXS[i]-u_g*scatXS[i])))
			

		#fills the transition matrix/scattering kernel
		for j in range(0,options.numGroups):
			scat[j,k-1]=M.data['fuel']['Ex'+str(j+1)][k]+M.data['moderator']['Ex'+str(j+1)][k]+M.data['poison']['Ex'+str(j+1)][k]
			#adding these moderator and poison values makes it unstable
			      
###############################################################################


	n=0
	while n<3:
		
		if n !=0:
			for i in range(0,options.numBins):
				totXS[i]=totXS[i]/2
				scatXS[i]=scatXS[i]/2
				diffcoef[i]=(1/(3*(totXS[i]-u_g*scatXS[i]))) 

		n=n+1
		
		#Builds the linear system of equations
		A=Construct()
		A.constructA(options, diffcoef, scat, totXS, numDensity)
		
		
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
