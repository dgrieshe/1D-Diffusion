#last edited by Miriam Rathbun on 8/15/2016
#This script solves the discrete diffusion equations



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np 
import matplotlib.pyplot as plt
from diffOpts import *
from plotter import * 
from GetFileName import *
from solver import *
from Construct import *
from nuclide import *
from material import *
from depletion import *

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
	nBins=options.numBins
	nGrps=options.numGroups
	
	#variables
	u_g=0.33
	n=0


	totXS=[]
	absXS=[]
	scatXS=[]
	fisXS=[]
	diffcoef=[]
	
	
	source=np.zeros(nBins*nGrps)
	#group-to-group scattering
	Gscat=np.zeros((nBins,nGrps*nGrps))
	sourceGroup=[]
	
	N=Nuclides()
	N.read()
	M=Material()
	M.read()
	M.updateNDarray(nBins,n)
	NDarray=M.NDarray
	
###############################################################
	
	while n<40:
		#first iteration: creates macroscopic cross section arrays
		#totXS, scatXS, and diffcoef with original number densities
		if n == 0:
			for k in range(1,nGrps+1):
				for i in range(0,nBins):
					totXS.append(N.data['fuel']['totxs'][k]*NDarray[i,0]+N.data['moderator']['totxs'][k]*NDarray[i,1]+N.data['poison']['totxs'][k]*NDarray[i,2])
					scatXS.append(N.data['fuel']['scatxs'][k]*NDarray[i,0]+N.data['moderator']['scatxs'][k]*NDarray[i,1]+N.data['poison']['scatxs'][k]*NDarray[i,2])
					fisXS.append(N.data['fuel']['fisxs'][k]*NDarray[i,0])
					diffcoef.append(1/(3*(totXS[i]-u_g*scatXS[i])))							

		else:
			D=Depletion()
			D.var()
			D.forEuler(sol.x*options.delta,NDarray)
			NDarray=D.NDarray
			
			plt.plot(NDarray[:,0])
			plt.savefig('./output/numdensity')
			
			#further iterations: updates totXS, scatXS, and diffcoef with new number densities
			#M.updateNDarray(nBins,n)
			#because the XS arrays have length nBins*nGrps, NDarray must iterate as i-nBins*(k-1)
			for k in range(1, nGrps+1):
				for i in range(nBins*(k-1),nBins*k):
					totXS[i]=N.data['fuel']['totxs'][k]*NDarray[i-nBins*(k-1),0]+N.data['moderator']['totxs'][k]*NDarray[i-nBins*(k-1),1]+N.data['poison']['totxs'][k]*NDarray[i-nBins*(k-1),2]
					scatXS[i]=N.data['fuel']['scatxs'][k]*NDarray[i-nBins*(k-1),0]+N.data['moderator']['scatxs'][k]*NDarray[i-nBins*(k-1),1]+N.data['poison']['scatxs'][k]*NDarray[i-nBins*(k-1),2]
					fisXS[i]=N.data['fuel']['fisxs'][k]*NDarray[i-nBins*(k-1),0]
					diffcoef[i]=(1/(3*(totXS[i]-u_g*scatXS[i])))
				
		#fills the group-to-group scattering/transition matrix per bin
		for i in range (0,nBins):
			count=1
			gcount=1
			for j in range(0,nGrps*nGrps):
				Gscat[i,j]=N.data['fuel']['Ex'+str(gcount)][count]*NDarray[i,0]+N.data['moderator']['Ex'+str(gcount)][count]*NDarray[i,1]+N.data['poison']['Ex'+str(gcount)][count]*NDarray[i,2]
				count=count+1
				if count == nGrps+1:
					count=1
					gcount=gcount+1
				
		
		#Builds the linear system of equations
		A=Construct()
		A.constructA(options, diffcoef, Gscat, totXS, NDarray)
		
		n=n+1
		
		#sets the initial value of the source
		for i in range (0,len(source)):
			source[i]=1/options.length

###############################################################################
#calculates the solution x to Ax=B 	


		A.invertA()
		sol=Solve()
		sol.solve(options, A.inv, source, fisXS, NDarray, N.data)
	
	
	results=Plotter()
	results.plot(sol.x,1,nBins,nGrps,name)
    
###############################################################################
