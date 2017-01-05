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
        self.energyPerFission = EperFission
        self.powerLevel = powerLevel
        # Cross sections
        self.fuel_absxs = N.data['fuel']['absxs'][1]
        self.fuel_fisxs = N.data['fuel']['fisxs'][1]
        self.nuclideList = N.nuclideList
        self.poisonList = N.poisonList
        # self.n is the number of columns in NDarray before 
            # the poisons start
        self.n = len(N.nuclideList)-len(N.poisonList)
        self.decayCST = []
        self.p_absxs = []
        self.Yield = []
        for p in self.poisonList:
            self.decayCST.append(N.data[p]['decayCST'])
            self.Yield.append(N.data[p]['yield'])
            self.p_absxs.append(N.data[p]['absxs'][1])
            
        
###############################################################
        
    def LocalEXP(self, flux, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("local matrixEXP")
        

        self.NDarray = NDarray

        ###################################
        # Step power normalization

        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))

        for i in range(0,len(NDarray)):
            summ = 0
            m = 0
            for p in self.poisonList:
                summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                m = m+1
            summation[i] = summ

        # Flux is already multiplied by delta. See input

        if PowerNormType == 'average':
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*self.powerLevel/sum(power)
            # Compute new power matrix to use in renormalization
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
        elif PowerNormType == 'explicit':
            # powerP1 is "part 1" of the power, where summation
                # is "part 2"
            powerP1[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]
            # I found that '+ summation[:]' and '- summation[:]'
                # do not cancel each other out
            #print sum(powerP1[:] + summation[:] - summation[:])
            #print sum(powerP1)
            flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
            # Compute new power matrix to use in renormalization
            power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
        #print sum(power)
        ###################################


        # Begin matrixEXP
        for i in range(0,len(NDarray)):
            #print("running...")

            # self.A is the matrix with the microscopic xs and 
                # decay constants
            # B is an identity matrix needed for the matrix
                # exponential calculation
            # self.NDbin is the matrix with the number densities
                # for a given bin i
            self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
            B = np.zeros((len(self.nuclideList),len(self.nuclideList)))
            self.NDbin = np.zeros(len(self.nuclideList))

            n = 0
            for nuclide in self.nuclideList:
                row = n
                self.NDbin[row] = self.NDarray[i,n]
                for col in range(0,len(self.A)):
                    if col == row:
                        self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                        B[row,col] = 1
                if n >= self.n:
                    self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.energyPerFission
                n = n+1



            # Depletion
            j = self.t
            while j <= self.t*self.num:

                if j == self.t:
                    self.NDarray[i,:] = np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))

                else:
                    # Power renormalization
                    # Flux is already multiplied by delta
                    # FUTURE WORK: other nuclides have fissxs
                    if PowerNormType == 'average':
                        poweri = flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                        flux[i] = flux[i]*power[i]/poweri
                    elif PowerNormType == 'explicit':
                        #print(power[i])
                        poweri = (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]+summation[i]
                        flux[i] = flux[i]*(power[i]-summation[i])/(poweri-summation[i])
                        #print (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*NumDensities[SSnum-1,0]+summation[i]
                    if flux[i] < 0:
                        print("Error: negative flux in bin %i" % i)



                    self.NDarray[i,:] = np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))
                    #for count in range(0,len(self.nuclideList)):
                    #    if np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))[count] < 0:
                    #        print("Error: negative number density in bin %i for nuclide %i" % (i,count))


                # Update self.NDbin between substeps
                self.NDbin[:] = np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))[:]

                # Update summation between substeps
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ

                j = j+self.t


            
###############################################################

    def GlobalEXP(self, flux, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("global matrixEXP")


        self.NDarray = NDarray
        #print("NDarray")
        #print self.NDarray

        # Step power normalization
        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))

        for i in range(0,len(NDarray)):
            summ = 0
            m = 0
            for p in self.poisonList:
                summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                m = m+1
            summation[i] = summ

        # Flux is already multiplied by delta.
        # I can use fisXS here because no depletion has yet 
            # occured

        if PowerNormType == 'average':
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*self.powerLevel/sum(power)
            # Compute new power matrix as a check
            #power[:] = flux[:]*self.energyPerFission*fisXS[:]
        elif PowerNormType == 'explicit':
            powerP1[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
            # Compute new power matrix as a check
            #power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
        #print sum(power)


        # Sub-stepping
        j = self.t
        while j <= self.t*self.num:
            #print("running")

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
                B = np.zeros((len(self.nuclideList),len(self.nuclideList)))
                self.NDbin = np.zeros(len(self.nuclideList))
    
                n = 0
                for nuclide in self.nuclideList:
                    row = n
                    self.NDbin[row] = self.NDarray[i,n]
                    for col in range(0,len(self.A)):
                        if col == row:
                            self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
                            B[row,col] = 1
                    if n >= self.n:
                        self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.energyPerFission
                    n = n+1

                # Depletion
                # Update number densities in bin i
                self.NDarray[i,:] = np.dot(self.NDbin,scipy.sparse.linalg.expm_multiply(self.A,B))
                # Update summation in bin i
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ
                # Update bin power
                # FUTURE WORK: not only the fuel has fisxs
                if PowerNormType == 'average':
                    power[i] = flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                elif PowerNormType == 'explicit':
                    powerP1[i] = (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]

            # Power renormalization
            if PowerNormType == 'average':
                flux[:] = flux[:]*self.powerLevel/sum(power)
                # Compute new power matrix as a check
                #power[:] = flux[:]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]
            elif PowerNormType == 'explicit':
                flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
                # Compute new power matrix as a check
                #power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]+summation[:]
            #print sum(power)


            # Next sub-step
            j = j+self.t
        
###############################################################

    def LocalEuler(self, flux, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("local forward Euler")
        print("Forward Euler is not accurate for long time steps!")

        self.NDarray = NDarray


        ###################################
        # Step power normalization

        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))

        for i in range(0,len(NDarray)):
            summ = 0
            m = 0
            for p in self.poisonList:
                summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                m = m+1
            summation[i] = summ


        # Flux is already multiplied by delta.

        if PowerNormType == 'average':
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*self.powerLevel/sum(power)
            # Compute new power matrix to use in renormalization
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
        elif PowerNormType == 'explicit':
            # powerP1 is "part 1" of the power, where summation
                # is "part 2"
            powerP1[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]
            # I found that '+ summation[:]' and '- summation[:]'
                # do not cancel each other out
            #print sum(powerP1[:] + summation[:] - summation[:])
            #print sum(powerP1)
            flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
            # Compute new power matrix to use in renormalization
            power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
        #print sum(power)
        ###################################



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
                    self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.energyPerFission
                n = n+1



            
            # Depletion
            j = self.t
            # There will be self.num number of sub-steps
            while j <= self.t*self.num:

                if j == self.t:
                    self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])

                else:
                    # Power renormalization
                    # Flux is already multiplied by delta
                    # See input for forEuler
                    # FUTURE WORK: realize that other
                        # nuclides could have a fission xs
                    if PowerNormType == 'average':
                        poweri = flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                        flux[i] = flux[i]*power[i]/poweri
                    elif PowerNormType == 'explicit':
                        #print(power[i])
                        poweri = (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]+summation[i]
                        flux[i] = flux[i]*(power[i]-summation[i])/(poweri-summation[i])
                        #print (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*NumDensities[SSnum-1,0]+summation[i]
                    if flux[i] < 0:
                        print("Error: negative flux in bin %i" % i)



                    self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])
                    for count in range(0,len(self.nuclideList)):
                        if np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][count] < 0:
                            print("Error: negative number density in bin %i for nuclide %i" % (i,count))


                # Update self.NDbin between substeps
                self.NDbin[:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][:]

                # Update summation between substeps
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ

                j = j+self.t



###############################################################

    def GlobalEuler(self, flux, NDarray, fisXS, YieldList, PowerNormType, N):
        #print("global forward Euler")
        print("Forward Euler is not accurate for long time steps!")


        self.NDarray = NDarray
        #print("NDarray")
        #print self.NDarray

        # Step power normalization
        power = np.zeros(len(NDarray))
        powerP1 = np.zeros(len(NDarray))
        summation = np.zeros(len(NDarray))

        for i in range(0,len(NDarray)):
            summ = 0
            m = 0
            for p in self.poisonList:
                summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                m = m+1
            summation[i] = summ

        # Flux is already multiplied by delta. See input 
            # for forEuler
        # I can use fisXS here because no depletion has yet 
            # occured

        if PowerNormType == 'average':
            power[:] = flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*self.powerLevel/sum(power)
            # Compute new power matrix as a check
            #power[:] = flux[:]*self.energyPerFission*fisXS[:]
        elif PowerNormType == 'explicit':
            powerP1[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]
            flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
            # Compute new power matrix as a check
            #power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
        #print sum(power)


        # Sub-stepping
        j = self.t
        while j <= self.t*self.num:

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
                        self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.energyPerFission
                    n = n+1

                # Depletion
                # Update number densities in bin i
                self.NDarray[i,:] = np.array([self.NDbin+j*np.dot(self.A,self.NDbin)])[0][:]
                # Update summation in bin i
                summ = 0
                m = 0
                for p in self.poisonList:
                    summ = summ + self.Yield[m]*self.decayCST[m]*self.NDarray[i,m+self.n]
                    m = m+1
                summation[i] = summ
                # Update bin power
                # FUTURE WORK: not only fuel has fisxs
                if PowerNormType == 'average':
                    power[i] = flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]
                elif PowerNormType == 'explicit':
                    powerP1[i] = (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[i,0]

            # Power renormalization
            if PowerNormType == 'average':
                flux[:] = flux[:]*self.powerLevel/sum(power)
                # Compute new power matrix as a check
                #power[:] = flux[:]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]
            elif PowerNormType == 'explicit':
                flux[:] = flux[:]*(self.powerLevel-sum(summation))/sum(powerP1)
                # Compute new power matrix as a check
                #power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*N.data['fuel']['fisxs'][1]*self.NDarray[:,0]+summation[:]
            #print sum(power)


            # Next sub-step
            j = j+self.t





###############################################################
# end class