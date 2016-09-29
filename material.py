# class containing material object. Replaces all the 
#microscopic cross section data from Nuclides with
#macroscopic cross section data by multiplying by 
#number densities read from an input file

import numpy as np
from nuclide import *

class Material():
	
###############################################################
	
	def __init__(self):
		self
###############################################################		
		
	def read(self):
		inpFile = open('Materials/NumDensities.inp','r')
		for line in inpFile:
			
			#Remove trailing whitespace
			line = line.strip()
			#Remove newline characters
			line = line.strip('\n')
			#Remove string after comment character (#)
			line, scratch1, scratch2 = line.partition('#')
			#Skip empty lines 
			if len(line) == 0:
				continue  
				
				
			if line == 'material' or line == 'end':
				continue
			else:
				keyword, arguments = line.split(' ',1)
				if keyword == 'fuel':
					self.NDfuel = float(arguments)
					
				elif keyword == 'moderator':
					self.NDmod = float(arguments)
					
				elif keyword == 'poison':
					self.NDpoison = float(arguments)

###############################################################
					
	def updateNDarray(self,nBins,n):
		if n==0:
			self.NDarray=np.zeros((nBins,3))
			for i in range(0,nBins):
				self.NDarray[i,0]=self.NDfuel
				self.NDarray[i,1]=self.NDmod
				self.NDarray[i,2]=self.NDpoison
		else:
			for i in range(0,nBins):
				self.NDarray[i,0]=self.NDarray[i,0]/2 #fuel
				self.NDarray[i,1]=self.NDarray[i,1]/2 #mod
				self.NDarray[i,2]=self.NDarray[i,2]/2 #poison
		
	
############################################################### 
#end
    	