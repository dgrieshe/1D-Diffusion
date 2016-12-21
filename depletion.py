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

	def var(self, N, powerLevel, nYield, EperFission):
		
		self.t = 0.0001
		self.num = 5
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
		#print self.n
		self.decayCST = []
		self.p_absxs = []
		self.Yield = []
		for p in self.poisonList:
			self.decayCST.append(N.data[p]['decayCST'])
			self.Yield.append(N.data[p]['yield'])
			self.p_absxs.append(N.data[p]['absxs'][1])
			
		
###############################################################
		
	def matrixEXP(self, flux, NDarray):
		
		# This method only accommodates 1 poison nuclide.
		# It was not updated for more than 1 poison.
		

		self.NDarray = NDarray
		for i in range(0,len(NDarray)):
			A = np.matrix([[-self.fuel_absxs*flux[i],0],[self.y*self.fuel_fisxs*flux[i],-(self.poison_absxs)*flux[i]-self.poison_decay]])
			B = np.matrix([[1,0],[0,1]])
			self.ND = np.matrix([NDarray[i,0],NDarray[i,2]])
			NumDensities=np.dot(self.ND,scipy.sparse.linalg.expm_multiply(A,B,start=0,stop=self.t*self.num,num=self.num))
			self.NDarray[i,0] = NumDensities[len(NumDensities)-1,0]
			self.NDarray[i,2] = NumDensities[len(NumDensities)-1,1] 
			
		
###############################################################

	def forEuler(self, flux, NDarray, fisXS, YieldList, PowerNormType, N):
		# FUTURE WORK: this part needs to be updated for 
		# boron burnup

		self.NDarray = NDarray
		#print("NDarray")
		#print self.NDarray

		# Step power normalization
		power = np.zeros(len(NDarray))
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
			# Compute new power matrix to use in renormalization
			power[:] = flux[:]*self.energyPerFission*fisXS[:]
		elif PowerNormType == 'explicit':
			power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
			flux[:] = flux[:]*(self.powerLevel-sum(summation))/(sum(power)-sum(summation))
			# Compute new power matrix to use in renormalization
			power[:] = (1-sum(YieldList))*flux[:]*self.energyPerFission*fisXS[:]+summation[:]
		#print sum(power)



        # Begin forEuler
		for i in range(0,len(NDarray)):

			# self.A is the matrix with the microscopic xs and 
			    # decay constants
			# self.ND is the matrix with the number densities
			self.A = np.zeros((len(self.nuclideList),len(self.nuclideList)))
			self.ND = np.zeros(len(self.nuclideList))

			n = 0
			for nuclide in self.nuclideList:
				row = n
				self.ND[row] = NDarray[i,n]
				for col in range(0,len(self.A)):
					if col == row:
						self.A[row,col] = -(N.data[nuclide]['absxs'][1]*flux[i]+N.data[nuclide]['decayCST'])
				if n >= self.n:
					self.A[row,0] = N.data[nuclide]['yield']*N.data['fuel']['fisxs'][1]*flux[i]*self.energyPerFission
				n = n+1


			
			# Sub-stepping
			j = self.t
			# SSnum is the sub-step number
			# There will be self.num number of sub-steps
			    # and SSnum keeps counts which one is occuring
			SSnum = 0

			while j <= self.t*self.num:


				if j == self.t:
					NumDensities = np.array([self.ND+j*np.dot(self.A,self.ND)])

				else:

					# Power renormalization
					# Flux is already multiplied by delta
					# See input for forEuler
					# I need to recalc fisXS because fuel 
					    # number density has changed
					if PowerNormType == 'average':
						poweri = flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*NumDensities[SSnum,0]
						flux[i] = flux[i]*power[i]/poweri
					elif PowerNormType == 'explicit':
						poweri = (1-sum(YieldList))*flux[i]*self.energyPerFission*N.data['fuel']['fisxs'][1]*NumDensities[SSnum-1,0]+summation[i]
						flux[i] = flux[i]*(power[i]-summation[i])/(poweri-summation[i])
					if flux[i] < 0:
						print("Error: negative flux in bin %i" % i)
						print(power[i]-summation[i])
						print(poweri-summation[i])



					NumDensities = np.concatenate((NumDensities,np.array([self.ND+j*np.dot(self.A,self.ND)])))
					for count in range(0,12):
					    if np.array([self.ND+j*np.dot(self.A,self.ND)])[0][count] < 0:
						    print("Error: negative number density in bin %i for nuclide %i" % (i,count))


				# Update self.ND between substeps
				self.ND[:] = np.array([self.ND+j*np.dot(self.A,self.ND)])[0][:]

				# Update summation between substeps
				summ = 0
				m = 0
				for p in self.poisonList:
					summ = summ + self.Yield[m]*self.decayCST[m]*NumDensities[SSnum,m+self.n]
					m = m+1
				summation[i] = summ

				SSnum = SSnum+1
				j = j+self.t


			# Updating NDarray
			# self.n is the number of columns in NDarray 
			    # before the poisons start
			self.NDarray[i,:] = NumDensities[self.num-1,:]

			#print self.NDarray[i,:]

		
###############################################################
# end class
