#last edited by Miriam Rathbun on 11/1/2016

class DiffusionOpts1D:
	
	def __init__(self):
		self.length=-1
		
	def read(self, filename):
		inpFile = open(filename, 'r')
		
		
		
		for line in inpFile:  


			#Remove trailing white space
			line = line.strip()
			#Remove newline characters
			line = line.strip('\n')
			#Remove string after comment character (#)
			line, scratch1, scratch2 = line.partition('#')
			#Skipe empty lines
			if len(line) == 0:
				continue


			keyword, arguments = line.split(' ', 1)
			if keyword == 'length':
				self.length = float(arguments)
				
			elif keyword == 'numgroups':
				self.numGroups = arguments
				self.numGroups = int(self.numGroups)
				
			elif keyword == 'numbins':
				self.numBins = arguments
				self.numBins = int(self.numBins)
				
			elif keyword == 'powerLevel':
				self.powerLevel = float(arguments)
				
			elif keyword == 'FneutronYield':
				self.nYield = arguments
				self.nYield = float(self.nYield)
				
			elif keyword == 'EperFission':
				self.EperFission = float(arguments)

			elif keyword == 'PowerNorm[average;explicit]':
				self.PowerNorm = arguments
				if (self.PowerNorm != 'average' and self.PowerNorm != 'explicit'):
					print('Error: PowerNorm entry is invalid. Use average or explicit. Delete all empty lines below it.')

			elif keyword == 'ConvergeError':
				self.ConvError = float(arguments)
				
			else:
				continue
					
						
		self.delta=self.length/self.numBins
							
								
                
                
				
				
				  