# Basic class containing an object for storing options for a simple 1D diffusion
# theory solver.

from errorHandler import *

class DiffusionOptions1D:
    # Initialization (constructor) routine.  This routine sets default values 
    # for all of the options.

    def __init__(self):
        self.length = -1

    # Routine for reading diffusion theory options from a file.

    def read(self, filename):
        inpFile = open(filename, 'r')

        for line in inpFile:
            # Remove any leading and trailing whitespace from the line.
            line = line.strip()

            # Eliminate any portion of the string that follows a comment
            # character (#)
            line, scratch1, scratch2 = line.partition('#')

            # If the line is empty the cycle to the next line.
            if len(line) == 0:
                continue

            # Perform casefolding on the line to eliminate any potential
            # differences in capitalization.
            line = line.casefold()

            # Split each line into a keyword and arguments string based on the
            # first space.
            keyword, arguments = line.split(' ', 1)
            
            if keyword == 'length':
                if(arguments.find(' ') > 0):
                    self.length, scratch = arguments.split(' ', 1)
                else:
                    self.length = arguments

                # Convert to a floating point value
                self.length = float(self.length)

                if(self.length <= 0):
                    fatalError("Invalid 1D slab length entered: "+str(self.length))

            elif keyword == 'numgroups':
                if(arguments.find(' ') > 0):
                    self.numGroups, scratch = arguments.split(' ', 1)
                else:
                    self.numGroups = arguments

                # Convert to a floating point value
                self.numGroups = int(self.numGroups)

                if(self.numGroups <= 0):
                    fatalError("Invalid number of energy groups: "+str(self.numGroups))

# Uncomment after support for all arguments has been added.
#            else:
#                  fatalError("Invalid input entry: "+keyword.strip())

        # Now we check to make sure that the necessary input parameters have 
        # all been read in.

        if self.length < 0:
          fatalError('1D slab length must be defined in the input file')
