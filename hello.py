import sys
#first time use: sudo apt-get install python-numpy 
import numpy as np
from diffusionOpts import *
from nuclide import *

print('Hello, World! (python)')
m=0
i=0
j=9
print(m)


while i==0:
	print(m)
	m=m+1
	if m==5:
		i=i+5

for j in range (0,3):
	print('hi')

print(np.cos(np.pi/6))
print(np.sqrt(3)/2)


###########################################################################
#
# Simple example.  Create a nuclide 'A' with name U235, then read in the
# xs data and print out the parameters.

A = Nuclide('U235')
A.read('input_filename')
print(A.name)
print(A.totalXs)

# Here is a slightly more complex example where we add a Nuclide object
# to a dictionary hashed by the nuclide name (e.g., the array key is the
# nuclide name.  Note that we can access all of the methods attached to
# the Nuclide object directly from the dictionary entry.

nuclides = {}
thisName = 'U238'
nuclides[thisName] = Nuclide(thisName)
nuclides[thisName].read('dummyFilename.txt')

print(nuclides['U238'].name)
print(nuclides['U238'].totalXs)

###########################################################################
#
# Example to read in diffusion options.


solverOptions = DiffusionOptions1D()
solverOptions.read('input1.inp')
print('Length:', solverOptions.length)
print('Groups:', solverOptions.numGroups)



