# class containing nuclide objects in a dictionary. 

class Nuclides:
    # Initialization (constructor) routine

    def __init__(self):
        self.name = 'name'

    # Routine for reading nuclide data from a file. 

    def read(self):
    	import numpy as np
    	inpFile = open('Materials/NuclideInput.inp', 'r')
    	self.data={} 
    	n=0
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
            
    		keyword, arguments = line.split(' ', 1)
    		if keyword == 'name':
    			self.name = arguments
    			
    		elif keyword == 'nuclide':
    			self.n = int(arguments)
    			
    		elif keyword == 'numGroups':
    			self.nGroups = int(arguments)
    			self.Gscat=np.zeros((self.nGroups,self.nGroups))
    			
    		elif keyword == 'totalXs':
    			self.totXs=[]
    			for i in range(1,self.nGroups+1):
    				self.totXs.append(float(line.split(' ',self.nGroups)[i]))
    				
    		elif keyword == 'absorptionXs':
    			self.absXs=[]
    			for i in range(1,self.nGroups+1):
    				self.absXs.append(float(line.split(' ',self.nGroups)[i]))
    				
    		elif keyword == 'fissionXs':
    			self.fisXs=[]
    			for i in range(1,self.nGroups+1):
    				self.fisXs.append(float(line.split(' ',self.nGroups)[i]))
    			
    		elif keyword == 'scatterXs':
    			self.scatXs=[]
    			for i in range(1,self.nGroups+1):
    				self.scatXs.append(float(line.split(' ',self.nGroups)[i]))
    			
    		elif keyword[:-1] == 'Gscat':  
    			for i in range(0, self.nGroups):
    				self.Gscat[n,i]=float(line.split(' ',self.nGroups)[i+1])
    			n=n+1
    
    			
    		if keyword == 'end':
    			n=0
    			self.d = {
    				self.name: {
    					'totXs' : {},
    					'absXs' : {},
    					'fisXs' : {},
    					'scatXs': {},
    				},
    			}
    			for i in range(1, self.nGroups+1):
    				self.d[self.name]['totXs'][i]=self.totXs[i-1]
    				self.d[self.name]['absXs'][i]=self.absXs[i-1]
    				self.d[self.name]['fisXs'][i]=self.fisXs[i-1]
    				self.d[self.name]['scatXs'][i]=self.scatXs[i-1]
    				self.d[self.name]['Ex'+str(i)]={}
    				for j in range(1,self.nGroups+1):
    					self.d[self.name]['Ex'+str(i)][j]=self.Gscat[i-1,j-1]
    			self.data.update(self.d)
    			
    			
    		else:
    			continue
    			

				

				
                
    	
        