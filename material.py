# class containing material object. Replaces all the 
#microscopic cross section data from Nuclides with
#macroscipic cross section data by multiplying by 
#number densities read from an input file

from nuclide import *

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
			self.data['fuel']['totXS'][i]=N.data['fuel']['totXS'][i]*self.NDfuel
			self.data['fuel']['absXS'][i]=N.data['fuel']['absXS'][i]*self.NDfuel
			self.data['fuel']['fisXS'][i]=N.data['fuel']['fisXS'][i]*self.NDfuel
			self.data['fuel']['scatXS'][i]=N.data['fuel']['scatXS'][i]*self.NDfuel
			self.data['moderator']['totXS'][i]=N.data['moderator']['totXS'][i]*self.NDmod
			self.data['moderator']['absXS'][i]=N.data['moderator']['absXS'][i]*self.NDmod
			self.data['moderator']['fisXS'][i]=N.data['moderator']['fisXS'][i]*self.NDmod
			self.data['moderator']['scatXS'][i]=N.data['moderator']['scatXS'][i]*self.NDmod
			self.data['poison']['totXS'][i]=N.data['poison']['totXS'][i]*self.NDpoison
			self.data['poison']['absXS'][i]=N.data['poison']['absXS'][i]*self.NDpoison
			self.data['poison']['fisXS'][i]=N.data['poison']['fisXS'][i]*self.NDpoison
			self.data['poison']['scatXS'][i]=N.data['poison']['scatXS'][i]*self.NDpoison
			for j in range(1,N.nGroups+1):
				self.data['fuel']['Ex'+str(i)][j]=N.data['fuel']['Ex'+str(i)][j]*self.NDfuel
				self.data['moderator']['Ex'+str(i)][j]=N.data['moderator']['Ex'+str(i)][j]*self.NDmod
				self.data['poison']['Ex'+str(i)][j]=N.data['poison']['Ex'+str(i)][j]*self.NDpoison

    			
    	