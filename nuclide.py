# Basic class containing a nuclide object and associated methods.

class Nuclide:
    # Initialization (constructor) routine

    def __init__(self, name):
        self.name = name

    # Routine for reading nuclide data from a file.  Note that these
    # cross section values are simply hard coded for now, but eventually
    # we will read them from the file specified on the 'filename' argument.

    def read(self):
    	inpFile = open('Materials/%s.inp' % self.name, 'r')
    	
    	for line in inpFile:  
			keyword, arguments = line.split(' ', 1)
			if keyword == 'name':
				self.name = arguments
				
			else:
				if keyword == 'totalXs':
					self.totalXs = arguments
					self.totalXs = float(self.totalXs)
					
				else:
					if keyword == 'absXs':
						self.absXs = arguments
						self.absXs = float(self.absXs)
						
					else:
						if keyword == 'fissXs':
							self.fissXs = arguments
							self.fissXs = float(self.fissXs)
							
						else:
							if keyword == 'scatXs':
								self.scatXs = arguments
								self.scatXs = float(self.scatXs)
								
							else:
								continue
                
    	
        