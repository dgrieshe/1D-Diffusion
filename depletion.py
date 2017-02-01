# class to compute depletion
# Last modified by Miriam Rathbun on 12-19-2016

import numpy as np
import scipy.sparse.linalg
from nuclide import *
from material import *

class Depletion():
    
###############################################################
    
    def __init__(self):
        self

###############################################################

    def var(self, N, t, num, powerLevel, nYield, EperFission):
        
        self.t = t
        self.num = num
        # Fission neutron yield
        self.y = nYield
        self.EperFission = EperFission
        self.powerLevel = powerLevel
        self.nuclideList = N.nuclideList
        self.poisonList = N.poisonList
        # self.n is the number of columns in NDarray before 
            # the poisons start
        self.n = len(N.nuclideList)-len(N.poisonList)
        self.decayCST = []
        self.Yield = []
        for p in self.poisonList:
            self.decayCST.append(N.data[p]['decayCST'])
            self.Yield.append(N.data[p]['yield'])
            
        
###############################################################
        
    def LocalEXP(self, flux, delta, summation, NDarray, fisXS, YieldList, PowerNormType, N, powerPlot, INDEX, TS):
        #print("local matrixEXP")
        
        self.NDarray = NDarray
        self.INDEX = INDEX

        ###################################

        # Setup for step power normalization
        power = np.zeros(len(NDarray))
        #summationSS = np.zeros(len(NDarray))

        if PowerNormType == 'average':
            power[:] = flux[:]*self.EperFission*fisXS[:]
        elif PowerNormType == 'explicit':
            power[:] = (1-sum(YieldList))*flux[:]*self.EperFission*fisXS[:]+summation[:]
        #print sum(power)

        ###################################
        #POWER = 0 

        # Begin matrixEXP
        for i in range(0,len(NDarray)):
            #print("running...")

            # self.A is the matrix with the microscopic xs and 
                # decay constants
            # self.NDbin is the matrix with the number densities
                # for a given bin i
            self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
            self.NDbin = np.zeros(len(self.nuclideList))

            n = 0
            for nuclide in self.nuclideList:
                row = n
                self.NDbin[row] = self.NDarray[i,n]
                for col in range(0,len(self.A)):
                    if col == row:
                        self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                if n >= self.n:
                    self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.EperFission
                n = n+1
            #if i == 320:
            #    print self.A


            ###################################
            # Depletion
            j = self.t/self.num
            SS = 1
            while j <= self.t:

                self.NDarray[i,:] = np.dot(scipy.sparse.linalg.expm(self.A*self.t/self.num),self.NDbin)

                if PowerNormType == 'average':
                    poweri = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                    flux[i] = flux[i]*power[i]/(poweri)
                    # Compute new power matrix as a check
                    #poweri = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                elif PowerNormType == 'explicit':
                    summ = 0
                    m = 0
                    for p in self.poisonList:
                        summ = summ + self.decayCST[m]*self.NDarray[i,m+self.n]
                        m = m+1
                    summation[i] = summ
                    powerP1i = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0] #+summationSS[i]
                    flux[i] = flux[i]*(power[i]-summation[i])/(powerP1i) #-summationSS[i])
                    # Compute new power matrix as a check
                    powerP1i = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                    #poweri = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]+summationSS[i]
                    if flux[i] < 0:
                         print("Error: negative flux in bin %i" % i)
                #if i == 320:
                #    print powerP1i+summation[i]
                #    print summation[i]
                #    print (flux[i])
                #    print self.NDarray[i,3]


                powerPlot[(TS-1)*(self.num+1)+SS] = powerPlot[(TS-1)*(self.num+1)+SS]+flux[i]*delta
                #if SS == 1:
                #    print powerPlot[(TS-1)*(self.num+1)+SS]


                # Update self.NDbin between substeps
                self.NDbin[:] = self.NDarray[i,:][:]


                j = j+self.t/self.num
                SS = SS+1

                ###################################
  
            #POWER = POWER + poweri
        #print POWER

        for i in range(0,self.num):
            self.INDEX = self.INDEX+1


            
###############################################################

    def GlobalEXP(self, flux, delta, NDarray, fisXS, YieldList, PowerNormType, N, powerPlot, INDEX):
        #print("global matrixEXP")


        self.NDarray = NDarray
        self.INDEX = INDEX
        #print("NDarray")
        #print self.NDarray

        # Used for step power re-normalization
        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))


        # Sub-stepping
        j = self.t/self.num
        while j <= self.t:
            #print("running")
            #print j

            for i in range (0,len(NDarray)):

                # self.A is the matrix with the microscopic xs 
                    # and decay constants
                # B is an identity matrix needed for the matrix
                    # exponential calculation
                # self.NDbin is the matrix with the number 
                    # densities for a given bin i
                # May be able to move this definition of 
                    #A and NDbin outside the loop
                self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
                self.NDbin = np.zeros(len(self.nuclideList))
    
                n = 0
                for nuclide in self.nuclideList:
                    row = n
                    self.NDbin[row] = self.NDarray[i,n]
                    for col in range(0,len(self.A)):
                        if col == row:
                            self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                    if n >= self.n:
                        self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.EperFission
                    n = n+1
                #if i == 320:
                #    print self.A

                # Depletion
                # Update number densities in bin i
                self.NDarray[i,:] = np.dot(scipy.sparse.linalg.expm(self.A*self.t/self.num),self.NDbin)
                #if i == 320:
                #    print self.NDarray[i,3]
                #self.NDarray[i,:] = np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))
                # Update summation in bin i
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ
                # Update bin power
                # FUTURE WORK: not only the fuel has fisxs
                if PowerNormType == 'average':
                    power[i] = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                elif PowerNormType == 'explicit':
                    powerP1[i] = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]

            # Power renormalization
            if PowerNormType == 'average':
                flux[:] = flux[:]*self.powerLevel/(sum(power)*delta)
                # Compute new power matrix as a check
                #power[:] = flux[:]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]
            elif PowerNormType == 'explicit':
                #flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
                flux[:] = flux[:]*(self.powerLevel-sum(summation)*delta)/(sum(powerP1)*delta)
                # Compute new power matrix as a check
                #power[:] = (1-sum(YieldList))*flux[:]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]+summation[:]
            #print sum(power)*delta

            #powerPlot.append(sum(flux)*delta)
            powerPlot[self.INDEX] = sum(flux)*delta
            self.INDEX = self.INDEX+1


            # Next sub-step
            j = j+self.t/self.num
        
###############################################################

    def LocalEuler(self, flux, delta, summation, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("local forward Euler")

        self.NDarray = NDarray

        power = np.zeros(len(NDarray))
        if PowerNormType == 'average':
            power[:] = flux[:]*self.EperFission*fisXS[:]
        elif PowerNormType == 'explicit':
            power[:] = (1-sum(YieldList))*flux[:]*self.EperFission*fisXS[:]+summation[:]
        #print sum(power)*delta


        # Begin forEuler
        for i in range(0,len(NDarray)):

            # self.A is the matrix with the microscopic xs and 
                # decay constants
            # self.NDbin is the matrix with the number densities
                # for a given bin i
            # May be able to move this definition of A and NDbin 
                # outside the loop
            self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
            self.NDbin = np.zeros(len(self.nuclideList))

            n = 0
            for nuclide in self.nuclideList:
                row = n
                self.NDbin[row] = self.NDarray[i,n]
                for col in range(0,len(self.A)):
                    if col == row:
                        self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                if n >= self.n:
                    self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.EperFission
                n = n+1



            
            # Depletion
            j = self.t/self.num
            # There will be self.num number of sub-steps
            while j <= self.t:

                if j == self.t:
                    self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])

                else:
                    # Power renormalization
                    # Flux is already multiplied by delta
                    # See input for forEuler
                    # FUTURE WORK: realize that other
                        # nuclides could have a fission xs
                    if PowerNormType == 'average':
                        poweri = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                        flux[i] = flux[i]*power[i]/(poweri)
                        # Compute new power matrix as a check
                        poweri = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                    elif PowerNormType == 'explicit':
                        #print(power[i])
                        poweri = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]+summation[i]
                        flux[i] = flux[i]*(power[i]-summation[i])/((poweri-summation[i]))
                        # Compute new power matrix as a check
                        poweri = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]+summation[i]
                    if flux[i] < 0:
                        print("Error: negative flux in bin %i" % i)
                    #if i == 320:
                    #    print poweri



                    self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])
                    #for count in range(0,len(self.nuclideList)):
                    #    if np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][count] < 0:
                    #        print("Error: negative number density in bin %i for nuclide %i" % (i,count))


                # Update self.NDbin between substeps
                self.NDbin[:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][:]

                # Update summation between substeps
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ

                j = j+self.t/self.num



###############################################################

    def GlobalEuler(self, flux, delta, summation, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("global forward Euler")


        self.NDarray = NDarray
        #print("NDarray")
        #print self.NDarray
        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))




        # Sub-stepping
        j = self.t/self.num
        while j <= self.t:



            for i in range (0,len(NDarray)):
                # self.A is the matrix with the microscopic xs 
                    # and decay constants
                # self.NDbin is the matrix with the number 
                    # densities for a given bin i
                # May be able to move this definition of 
                    #A and NDbin outside the loop
                self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
                self.NDbin = np.zeros(len(self.nuclideList))
    
                n = 0
                for nuclide in self.nuclideList:
                    row = n
                    self.NDbin[row] = self.NDarray[i,n]
                    for col in range(0,len(self.A)):
                        if col == row:
                            self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                    if n >= self.n:
                        self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.EperFission
                    n = n+1



                # Depletion
                # Update number densities in bin i
                self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][:]
                # Update summation in bin i
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ
                # Update bin power
                # FUTURE WORK: not only fuel has fisxs
                if PowerNormType == 'average':
                    power[i] = flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                elif PowerNormType == 'explicit':
                    powerP1[i] = (1-sum(YieldList))*flux[i]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]



            # Power renormalization
            if PowerNormType == 'average':
                flux[:] = flux[:]*self.powerLevel/(sum(power)*delta)
                # Compute new power matrix as a check
                #power[:] = flux[:]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]
            elif PowerNormType == 'explicit':
                flux[:] = flux[:]*(self.powerLevel-sum(summation)*delta)/(sum(powerP1)*delta)
                # Compute new power matrix as a check
                #power[:] = (1-sum(YieldList))*flux[:]*self.EperFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]+summation[:]
            #print sum(power)*delta



            # Next sub-step
            j = j+self.t/self.num



###############################################################
# end class