#last edited by Miriam Rathbun on 7/9/2016

class DiffusionOpts1D:
	
	def __init__(self):
		self.length=-1
		
	def read(self, filename):
		inpFile = open(filename, 'r')
		
		
		
		for line in inpFile:  
			keyword, arguments = line.split(' ', 1)
			if keyword == 'length':
				self.length = arguments
				self.length = float(self.length)
				
			elif keyword == 'numgroups':
				self.numGroups = arguments
				self.numGroups = int(self.numGroups)
				
			elif keyword == 'numbins':
				self.numBins = arguments
				self.numBins = int(self.numBins)
				
			else:
				continue
					
						
		self.delta=self.length/self.numBins
							
								
                
                
				
				
				  