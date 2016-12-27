# Last edited by Miriam Rathbun on 12/19/2016
# This script solves the discrete diffusion equations



# Imports
# Use sudo apt-get install python.numpy if needed
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
    # Retains only the number associated with the input
        # file: "input1" becomes "1" This makes the 
        # output name nice
    name = f[len(f)-5]
    
    # read input
    options = DiffusionOpts1D()
    options.read(f)
    nBins = options.numBins
    nGrps = options.numGroups

    
    # Variables
    n = 0

    totXS = []
    absXS = []
    scatXS = []
    fisXS = []
    diffcoef = []
    
    
    source = np.zeros(nBins*nGrps)
    # Group-to-group scattering
    Gscat = np.zeros((nBins,nGrps*nGrps))
    sourceGroup = []
    #print Gscat
    
    N = Nuclides()
    N.read()
    M = Material()
    M.read()
    M.createNDarray(nBins,n)
    NDarray = M.NDarray
    #print NDarray
    #print len(N.nuclideList)-len(N.poisonList)
    #print N.YieldList
    
    
###############################################################
    
    while n<10:
        # First iteration: creates macroscopic cross section
            # arrays totXS, scatXS, and diffcoef with original
            # number densities

        # FUTURE WORK: make subroutine that updates totXS,
            # scatXS, fisXS, and diffcoef
        # FUTURE WORK: collapse the n=0 and n!=0 loops
        if n == 0:
            for k in range(1,nGrps+1):
                for i in range(0,nBins):
                    tot = 0
                    scat = 0
                    fis = 0
                    ugTOP = 0
                    ugBOT = 0
                    j = 0
                    for nuclide in N.nuclideList:
                        tot = tot + N.data[nuclide]['totxs'][k]*NDarray[i,j]
                        scat = scat + N.data[nuclide]['scatxs'][k]*NDarray[i,j]
                        fis = fis + N.data[nuclide]['fisxs'][k]*NDarray[i,j]
                        ugTOP = ugTOP + N.data[nuclide]['ug']*N.data[nuclide]['scatxs'][k]*NDarray[i,j]
                        ugBOT = ugBOT + N.data[nuclide]['scatxs'][k]*NDarray[i,j]
                        j = j+1
                    totXS.append(tot)
                    scatXS.append(scat)
                    fisXS.append(fis)
                    u_g = ugTOP/ugBOT
                    diffcoef.append(1/(3*(totXS[i]-u_g*scatXS[i])))



        else:
            D=Depletion()
            D.var(N, options.powerLevel, options.nYield, options.EperFission)
            if options.DepletionType == 'local':
                D.LocalEuler(sol.x*options.delta, NDarray, fisXS, N.YieldList, options.PowerNorm, N)
            elif options.DepletionType == 'global':
                D.GlobalEuler(sol.x*options.delta, NDarray, fisXS, N.YieldList, options.PowerNorm, N)
            NDarray = D.NDarray
            #print NDarray
            
            
            plt.plot(NDarray[:,0])
            plt.savefig('./output/numdensity')
            
            total = totXS
            # Further iterations: updates totXS, scatXS, and
                # diffcoef with new number densities Because the 
                # XS arrays have length nBins*nGrps, NDarray must 
                # iterate as i-nBins*(k-1)
            for k in range(1, nGrps+1):
                for i in range(nBins*(k-1),nBins*k):

                    # Reset local variables
                    tot = 0 
                    scat = 0
                    fis = 0
                    ugTOP = 0
                    ugBOT = 0
                    j = 0
                    for nuclide in N.nuclideList:
                        tot = tot + N.data[nuclide]['totxs'][k]*NDarray[i-nBins*(k-1),j]
                        scat = scat + N.data[nuclide]['scatxs'][k]*NDarray[i-nBins*(k-1),j]
                        fis = fis + N.data[nuclide]['fisxs'][k]*NDarray[i-nBins*(k-1),j]
                        ugTOP = ugTOP + N.data[nuclide]['ug']*N.data[nuclide]['scatxs'][k]*NDarray[i-nBins*(k-1),j]
                        ugBOT = ugBOT + N.data[nuclide]['scatxs'][k]*NDarray[i-nBins*(k-1),j]
                        j = j+1
                    totXS[i] = tot
                    scatXS[i] = scat
                    fisXS[i] = fis
                    u_g = ugTOP/ugBOT
                    diffcoef[i] = (1/(3*(totXS[i]-u_g*scatXS[i])))

                
        # Fills/Updates the group-to-group scattering/transition 
            # matrix per bin
        for i in range(0,nBins):
            count = 1
            gcount = 1
            for g in range(0,nGrps*nGrps):
                j = 0
                GscatXS = 0
                for nuclide in N.nuclideList:
                    GscatXS = GscatXS + N.data[nuclide]['Ex'+str(gcount)][count]*NDarray[i,j]
                    #print N.data[nuclide]['Ex'+str(gcount)][count]
                    #print NDarray[i,j]
                    #print GscatXS
                    j = j+1
                Gscat[i,g] = GscatXS
                count = count+1
                if count == nGrps+1:
                    count = 1
                    gcount = gcount+1
        #print Gscat
                
        
        # Builds the linear system of equations
        A=Construct()
        A.constructA(options, diffcoef, Gscat, totXS, NDarray)
        
        n = n+1
        
        # Sets the initial value of the source
        source[:] = 1/options.length

####################################################################
# Calculates the solution x to Ax=B     


        A.invertA()
        sol = Solve()
        sol.solve(options, A.inv, source, fisXS, NDarray, N.data, n)
    results = Plotter()
    results.plot(sol.x,1,nBins,nGrps,name,n)
    
###################################################################
