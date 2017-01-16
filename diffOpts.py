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
                self.numGroups = int(arguments)
                
            elif keyword == 'numbins':
                self.numBins = int(arguments)
                
            elif keyword == 'powerLevel':
                self.powerLevel = float(arguments)
                
            elif keyword == 'FneutronYield':
                self.nYield = float(arguments)
                
            elif keyword == 'EperFission':
                self.EperFission = float(arguments)

            elif keyword == 'PowerNorm[average;explicit]':
                self.PowerNorm = arguments
                if (self.PowerNorm != 'average' and self.PowerNorm != 'explicit'):
                    print('Error: PowerNorm entry is invalid. Use average or explicit.')
                    
            elif keyword == 'RenormType[local;global]':
                self.RenormType = arguments
                if (self.RenormType != 'local' and self.RenormType != 'global'):
                    print('Error: RenormType entry is invalid. Use local or global.')

            elif keyword == 'DepletionType[matrixEXP;forEuler]':
                self.DepletionType = arguments
                if (self.DepletionType != 'matrixEXP' and self.DepletionType != 'forEuler'):
                    print('Error: DepletionType entry is invalid. Use matrixEXP or forEuler.')

            elif keyword == 'ConvergeError':
                self.ConvError = float(arguments)

            elif keyword == 'numTimeStep':
                self.n = int(arguments)

            elif keyword == 'timeStep[sec]':
                self.timeStep = float(arguments)

            elif keyword == 'numSubStep':
                self.numSubStep = int(arguments)
                
            else:
                continue
                    
                        
        self.delta=self.length/self.numBins
                            
                                
                
                
                
                
                  