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
	
	#read input
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
	fuel_totXS=[]
	fuel_scatXS=[]
	fuel_DC=[]
	mod_totXS=[]
	mod_scatXS=[]
	mod_DC=[]
	pois_totXS=[]
	pois_scatXS=[]
	pois_DC=[]
	
	
	numDensity=np.zeros((options.numBins,3))
	numDensity[:,0]=M.NDfuel
	numDensity[:,1]=M.NDmod
	numDensity[:,2]=M.NDpoison
	source=np.zeros(options.numBins*options.numGroups)
	scat=np.zeros((options.numGroups,options.numGroups))
	sourceGroup=[]
	
	
###############################################################################	
	
#create slab that's all U-235
#loop through the groups
#fill other items which depend on groups, see below 

    

	for k in range(1,options.numGroups+1):
		for i in range(0,options.numBins):
			totXS.append(M.data['fuel']['totXs'][k]+M.data['moderator']['totXs'][k]+M.data['poison']['totXs'][k])
			scatXS.append(M.data['fuel']['scatXs'][k]+M.data['moderator']['scatXs'][k]+M.data['poison']['scatXs'][k])
			absXS.append(M.data['fuel']['absXs'][k]+M.data['moderator']['absXs'][k]+M.data['poison']['absXs'][k])
			fisXS.append(M.data['fuel']['fisXs'][k]+M.data['moderator']['fisXs'][k]+M.data['poison']['fisXs'][k])
			fuel_totXS.append(M.data['fuel']['totXs'][k])
			fuel_scatXS.append(M.data['fuel']['scatXs'][k])
			
			diffcoef.append(1/(3*(fuel_totXS[i]-u_g*fuel_scatXS[i])))
			
		#fills the transition matrix/scattering kernel
		for j in range(0,options.numGroups):
			scat[j,k-1]=M.data['fuel']['Ex'+str(j+1)][k]	
			      
###############################################################################


	n=0
	while n<3:
		
		if n !=0:
			for i in range(0,options.numBins):
				fuel_totXS[i]=fuel_totXS[i]/numDensity[i,0]
				fuel_scatXS[i]=fuel_scatXS[i]/numDensity[i,0]
			for i in range(0,options.numBins):
				for j in range(0,3):
					numDensity[i,j]=numDensity[i,j]/2
			for i in range(0,options.numBins):
				fuel_totXS[i]=fuel_totXS[i]*numDensity[i,0]
				fuel_scatXS[i]=fuel_scatXS[i]*numDensity[i,0]
				diffcoef[i]=(1/(3*(fuel_totXS[i]-u_g*fuel_scatXS[i]))) 
#the scattering kernel gets multiplied by numDensity in constructA
#for k in range(1,options.numGroups+1):
#	for j in range(0,options.numGroups):
#		scat[j,k-1]=scat[j,k-1]/2
		n=n+1
		
		#Builds the linear system of equations
		A=Construct()
		A.constructA(options, diffcoef, scat, fuel_totXS, numDensity)
		
		
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
