# Basic class containing a nuclide object and associated methods.

class Nuclide:
    # Initialization (constructor) routine

    def __init__(self, name):
        self.name = name

    # Routine for reading nuclide date from a file.  Note that these
    # cross section values are simply hard coded for now, but eventually
    # we will read them from the file specified on the 'filename' argument.

    def read(self, filename):
        self.totalXs = 10
        self.absXs   = 5
        self.fissXs  = 3
        self.scatXs  = 5
