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
    			self.totxs=[]
    			for i in range(1,self.nGroups+1):
    				self.totxs.append(float(line.split(' ',self.nGroups)[i]))
    				
    		elif keyword == 'absorptionXs':
    			self.absxs=[]
    			for i in range(1,self.nGroups+1):
    				self.absxs.append(float(line.split(' ',self.nGroups)[i]))
    				
    		elif keyword == 'fissionXs':
    			self.fisxs=[]
    			for i in range(1,self.nGroups+1):
    				self.fisxs.append(float(line.split(' ',self.nGroups)[i]))
    			
    		elif keyword == 'scatterXs':
    			self.scatxs=[]
    			for i in range(1,self.nGroups+1):
    				self.scatxs.append(float(line.split(' ',self.nGroups)[i]))
    			
    		elif keyword[:-1] == 'Gscat':  
    			for i in range(0, self.nGroups):
    				self.Gscat[n,i]=float(line.split(' ',self.nGroups)[i+1])
    			n=n+1
    
    			
    		if keyword == 'end':
    			n=0
    			self.d = {
    				self.name: {
    					'totxs' : {},
    					'absxs' : {},
    					'fisxs' : {},
    					'scatxs': {},
    				},
    			}
    			for i in range(1, self.nGroups+1):
    				self.d[self.name]['totxs'][i]=self.totxs[i-1]
    				self.d[self.name]['absxs'][i]=self.absxs[i-1]
    				self.d[self.name]['fisxs'][i]=self.fisxs[i-1]
    				self.d[self.name]['scatxs'][i]=self.scatxs[i-1]
    				self.d[self.name]['Ex'+str(i)]={}
    				for j in range(1,self.nGroups+1):
    					self.d[self.name]['Ex'+str(i)][j]=self.Gscat[i-1,j-1]
    			self.data.update(self.d)
    			
    			
    		else:
    			continue
    			

				

				
                
    	
        