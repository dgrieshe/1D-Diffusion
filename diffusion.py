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

# Variables for plotting
flux1 = []
flux2 = []
inp = 1

for f in fn.listfn:
    print f

    # read input
    options = DiffusionOpts1D()
    options.read(f)
    nBins = options.numBins
    nGrps = options.numGroups

    
    # Variables
    n = 0
    time = 0

    totXS = []
    absXS = []
    scatXS = []
    fisXS = []
    diffcoef = []
    massU = []
    timeMU = []
    
    
    source = np.zeros(nBins*nGrps)
    # Group-to-group scattering
    Gscat = np.zeros((nBins,nGrps*nGrps))
    sourceGroup = []
    powerPlot = np.zeros(1+(options.n-1)*(options.numSubStep+1))
    timePPlot = []
    pPlotINDEX = 0
    #print 1+(options.n-1)*(options.numSubStep+1)
    
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
    
    while n<options.n:
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
            massU.append(sum(NDarray[:,0])*options.delta)
            timeMU.append(time)
            massU1 = sum(NDarray[:,0])*options.delta
            #print massU[n]



        else:
            D=Depletion()
            D.var(N, options.timeStep, options.numSubStep, options.powerLevel, options.nYield, options.EperFission, time, timePPlot)
            if options.RenormType == 'local':
                if options.DepletionType == 'matrixEXP':
                    D.LocalEXP(flux, options.delta, summation, NDarray, fisXS, N.YieldList, options.PowerNorm, N, powerPlot, pPlotINDEX, massU, timeMU, n)
                    pPlotINDEX = D.INDEX
                elif options.DepletionType == 'forEuler':
                    D.LocalEuler(flux, options.delta, summation, NDarray, fisXS, N.YieldList, options.PowerNorm, N)
            elif options.RenormType == 'global':
                if options.DepletionType == 'matrixEXP':
                    D.GlobalEXP(flux, options.delta, NDarray, fisXS, N.YieldList, options.PowerNorm, N, powerPlot, pPlotINDEX, massU, timeMU)
                    pPlotINDEX = D.INDEX
                elif options.DepletionType == 'forEuler':
                    D.GlobalEuler(flux, options.delta, summation, NDarray, fisXS, N.YieldList, options.PowerNorm, N)
            NDarray = D.NDarray
            time = D.time
            timePPlot = D.timePPlot
            massU = D.massU
            timeMU = D.timeMU
            #print NDarray[320][0]
            massU.append(sum(NDarray[:,0])*options.delta)
            timeMU.append(time)
            #print massU[n]
            
            
            #plt.plot(NDarray[:,0])
            #plt.savefig('./output/numdensity')
            
            total = totXS
            # Further iterations: updates totXS, scatXS, and
                # diffcoef with new number densities because the 
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
        #print("solving1")
        sol = Solve()
        sol.solve(options, A.inv, source, fisXS, NDarray, N.data, n)
        #print("solving2")


####################################################################
# Prepares plotting

        # Power normalization to get "true" flux to plot
        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))
        flux = sol.x

        nn = len(N.nuclideList)-len(N.poisonList)
        decayCST = []
        p_absxs = []
        Yield = []
        for p in N.poisonList:
            decayCST.append(N.data[p]['decayCST'])
            Yield.append(N.data[p]['yield'])
            p_absxs.append(N.data[p]['absxs'][1])



        if n == 0:
            for k in range(1,nGrps+1):
                for i in range(0,nBins):
                    fis = 0
                    j = 0
                    for nuclide in N.nuclideList:
                        fis = fis + N.data[nuclide]['fisxs'][k]*NDarray[i,j]
                        j = j+1
                    fisXS.append(fis)

        for i in range(0,len(NDarray)):
            summ = 0
            m = 0
            for p in N.poisonList:
                summ = summ + decayCST[m]*NDarray[i,m+nn]
                m = m+1
            summation[i] = summ


        if options.PowerNorm == 'average':
            power[:] = flux[:]*options.EperFission*fisXS[:]
            flux[:] = flux[:]*options.powerLevel/(sum(power)*options.delta)
            # Compute new power matrix as a check
            #power[:] = flux[:]*options.EperFission*fisXS[:]
        elif options.PowerNorm == 'explicit':
            powerP1[:] = (1-sum(N.YieldList))*flux[:]*options.EperFission*fisXS[:]
            flux[:] = flux[:]*(options.powerLevel-sum(summation)*options.delta)/(sum(powerP1)*options.delta)
            # Compute new power matrix as a check
            #power[:] = (1-sum(N.YieldList))*flux[:]*options.EperFission*fisXS[:]+summation[:]
        #print sum(power)*options.delta

        #Storing fluxes for FluxRMS calc
        if inp == 1:
            flux1.append(sol.x)
        if inp == 2:
            flux2.append(sol.x)

        #powerPlot.append(sum(flux)*options.delta)
        powerPlot[pPlotINDEX] = sum(flux)*options.delta
        timePPlot.append(time)
        #print pPlotINDEX
        #print powerPlot[pPlotINDEX]
        pPlotINDEX = pPlotINDEX+1

####################################################################
# Plotting 

    results = Plotter()
    #results.plotFLUX(sol.x,1,nBins,nGrps,inp,n)
    #results.plotINTflux(timePPlot, powerPlot, inp)
    results.plotMASSU(timeMU, massU, massU1, n)
    inp = inp+1

    
#if inp == 2:
    #results.plotFluxRMS(sol.x, inp, n, flux1, flux2)
    
###################################################################
