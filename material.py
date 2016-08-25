# class containing material object. Replaces all the 
#microscopic cross section data from Nuclides with
#macroscipic cross section data by multiplying by 
#number densities read from an input file

from new_nuclide import *

class Material:
	
	def __init__(self):
		self
		
		
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
					
	def calc_macro(self):
		N=Nuclides()
		N.read()
		self.data=N.data
		for i in range(1, N.nGroups+1):
			self.data['fuel']['totXs'][i]=N.data['fuel']['totXs'][i]*self.NDfuel
			self.data['fuel']['absXs'][i]=N.data['fuel']['absXs'][i]*self.NDfuel
			self.data['fuel']['fisXs'][i]=N.data['fuel']['fisXs'][i]*self.NDfuel
			self.data['fuel']['scatXs'][i]=N.data['fuel']['scatXs'][i]*self.NDfuel
			self.data['moderator']['totXs'][i]=N.data['fuel']['totXs'][i]*self.NDmod
			self.data['moderator']['absXs'][i]=N.data['fuel']['absXs'][i]*self.NDmod
			self.data['moderator']['fisXs'][i]=N.data['fuel']['fisXs'][i]*self.NDmod
			self.data['moderator']['scatXs'][i]=N.data['fuel']['scatXs'][i]*self.NDmod
			self.data['poison']['totXs'][i]=N.data['fuel']['totXs'][i]*self.NDpoison
			self.data['poison']['absXs'][i]=N.data['fuel']['absXs'][i]*self.NDpoison
			self.data['poison']['fisXs'][i]=N.data['fuel']['fisXs'][i]*self.NDpoison
			self.data['poison']['scatXs'][i]=N.data['fuel']['scatXs'][i]*self.NDpoison
			for j in range(1,N.nGroups+1):
    					self.data['fuel']['Ex'+str(i)][j]=N.data['fuel']['Ex'+str(i)][j]*self.NDfuel
    					self.data['moderator']['Ex'+str(i)][j]=N.data['moderator']['Ex'+str(i)][j]*self.NDmod
    					self.data['poison']['Ex'+str(i)][j]=N.data['poison']['Ex'+str(i)][j]*self.NDpoison
    		
    			
    	