# class to compute depletion

import numpy as np
import scipy.sparse.linalg
from nuclide import *
from material import *

class Depletion():
	
###############################################################
	
	def __init__(self):
		self

###############################################################

	def var(self):
		
		N=Nuclides()
		N.read()
		
		self.t=0.001
		self.num=5
		#fission neutron yield
		self.y=2.3
		#cross sections
		self.fuel_absxs=N.data['fuel']['absxs'][1]
		self.fuel_fisxs=N.data['fuel']['fisxs'][1]
		self.poison_absxs=N.data['poison']['absxs'][1]
		#flux
		self.flux=10
		
###############################################################
		
	def matrixEXP(self,flux,NDarray):
		

		self.NDarray=NDarray
		for i in range(0,len(NDarray)):
			A=np.matrix([[-self.fuel_absxs*flux[i,0],0],[self.y*self.fuel_fisxs*flux[i,0],-(self.poison_absxs)*flux[i,0]]])
			B=np.matrix([[1,0],[0,1]])
			self.ND=np.matrix([NDarray[i,0],NDarray[i,2]])
			NumDensities=np.dot(self.ND,scipy.sparse.linalg.expm_multiply(A,B,start=0,stop=self.t*self.num,num=self.num))
			self.NDarray[i,0]=NumDensities[len(NumDensities)-1,0]
			self.NDarray[i,2]=NumDensities[len(NumDensities)-1,1] 
			
		
###############################################################

	def forEuler(self, flux, NDarray):
		
		self.NDarray=NDarray                                    
		for i in range(0,len(NDarray)):
			A=np.matrix([[-self.fuel_absxs*flux[i,0],0],[self.y*self.fuel_fisxs*flux[i,0],-(self.poison_absxs)*flux[i,0]]])
			self.ND=np.matrix([NDarray[i,0],NDarray[i,2]])
			NumDensities=np.array(self.ND+self.t*np.dot(self.ND,A))
			self.NDarray[i,0]=NumDensities[0,0]
			self.NDarray[i,2]=NumDensities[0,1] 
		
		
	#	i=self.t                                  
	#	while i<=self.t*self.num:                     
	#		if i==self.t:                                      
	#			self.NumDensities=np.array(self.ND+i*np.dot(self.ND,A))
	#		else:
	#			self.NumDensities=np.concatenate((self.NumDensities,np.array(self.ND+i*np.dot(self.ND,A))))
	#		i=i+self.t
		
###############################################################