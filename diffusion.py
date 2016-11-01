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

fn = FileName()
fn.GetFileName()

for f in fn.listfn:
	print f
	#retains only the number associated with the input file. Ex: "input1" becomes "1" 
	#This was made to make the output name nice
	name = f[len(f)-5]
	
	#read input
	options = DiffusionOpts1D()
	options.read(f)
	nBins = options.numBins
	nGrps = options.numGroups
	
	#variables
	u_g = 0.33
	n = 0


	totXS = []
	absXS = []
	scatXS = []
	fisXS = []
	diffcoef = []
	
	
	source = np.zeros(nBins*nGrps)
	#group-to-group scattering
	Gscat = np.zeros((nBins,nGrps*nGrps))
	sourceGroup = []
	
	N = Nuclides()
	N.read()
	M = Material()
	M.read()
	M.createNDarray(nBins,n)
	NDarray = M.NDarray
	#print NDarray
	#print N.nuclideList
	
###############################################################
	
	while n<15:
		#first iteration: creates macroscopic cross section arrays
		#totXS, scatXS, and diffcoef with original number densities
		if n == 0:
			for k in range(1,nGrps+1):
				for i in range(0,nBins):
					tot = 0
					scat = 0
					j = 0
					for nuclide in N.nuclideList:
						tot=tot+N.data[nuclide]['totxs'][k]*NDarray[i,j]
						scat=scat+N.data[nuclide]['scatxs'][k]*NDarray[i,j]
						j = j+1
					totXS.append(tot)
					scatXS.append(scat)
					#totXS.append(N.data['fuel']['totxs'][k]*NDarray[i,0]+N.data['moderator']['totxs'][k]*NDarray[i,1]+N.data['poison']['totxs'][k]*NDarray[i,2])
					#scatXS.append(N.data['fuel']['scatxs'][k]*NDarray[i,0]+N.data['moderator']['scatxs'][k]*NDarray[i,1]+N.data['poison']['scatxs'][k]*NDarray[i,2])
					fisXS.append(N.data['fuel']['fisxs'][k]*NDarray[i,0])
					diffcoef.append(1/(3*(totXS[i]-u_g*scatXS[i])))	

		else:
			D=Depletion()
			D.var(N, options.powerLevel, options.nYield, options.EperFission)
			D.forEuler(sol.x*options.delta, NDarray, fisXS)
			NDarray = D.NDarray
			#print NDarray
			
			plt.plot(NDarray[:,0])
			plt.savefig('./output/numdensity')
			
			total = totXS
			#further iterations: updates totXS, scatXS, and diffcoef with new number densities
			#because the XS arrays have length nBins*nGrps, NDarray must iterate as i-nBins*(k-1)
			for k in range(1, nGrps+1):
				for i in range(nBins*(k-1),nBins*k):
					tot = 0 
					scat = 0
					j = 0
					for nuclide in N.nuclideList:
						tot = tot+N.data[nuclide]['totxs'][k]*NDarray[i-nBins*(k-1),j]
						scat = scat+N.data[nuclide]['scatxs'][k]*NDarray[i-nBins*(k-1),j]
						j = j+1
					totXS[i] = tot
					scatXS[i] = scat
					#totXS[i]=N.data['fuel']['totxs'][k]*NDarray[i-nBins*(k-1),0]+N.data['moderator']['totxs'][k]*NDarray[i-nBins*(k-1),1]+N.data['poison']['totxs'][k]*NDarray[i-nBins*(k-1),2]
					#scatXS[i]=N.data['fuel']['scatxs'][k]*NDarray[i-nBins*(k-1),0]+N.data['moderator']['scatxs'][k]*NDarray[i-nBins*(k-1),1]+N.data['poison']['scatxs'][k]*NDarray[i-nBins*(k-1),2]
					fisXS[i] = N.data['fuel']['fisxs'][k]*NDarray[i-nBins*(k-1),0]
					diffcoef[i] = (1/(3*(totXS[i]-u_g*scatXS[i])))
				
		#fills the group-to-group scattering/transition matrix per bin
		for i in range (0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				j = 0
				for nuclides in N.nuclideList:
					Gscat[i,g] = N.data[nuclide]['Ex'+str(gcount)][count]*NDarray[i,j]
					j = j+1
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1
				
		
		#Builds the linear system of equations
		A=Construct()
		A.constructA(options, diffcoef, Gscat, totXS, NDarray)
		
		n = n+1
		
		#sets the initial value of the source
		for i in range (0,len(source)):
			source[i] = 1/options.length

###############################################################################
#calculates the solution x to Ax=B 	


		A.invertA()
		sol = Solve()
		sol.solve(options, A.inv, source, fisXS, NDarray, N.data)
	
	
	results = Plotter()
	results.plot(sol.x,1,nBins,nGrps,name)
    
###############################################################################
