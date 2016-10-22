#class created to solve Ax=B

class Solve:
	
	def __init__(self):
		self
		
	def solve(self, options, Ainv, source, fisXS, NDarray, Ndata):
		
		import numpy as np

		# Create array variables.
		lastSource = np.zeros(options.numBins*options.numGroups)
		errorDiff = np.zeros(options.numBins*options.numGroups)

		# Initialize local variables
		error=1000 
		j=0
                
                # Loop until the source difference between iterations falls
                # below the specified threshold.
                
                # Note: we should probably read the threshold from the input
                # file
                
		while error > 10 ** (-5):
                	# Reset local variables
			k = 0
			j = j+1
                        errorDiff[:] = 0.
                        
                        # Save the previous source so that we can compute the
                        # residual error for the iteration.
                        
                        lastSource[:] = source[:]
                        
			# Take the dot product of A-inverse and B to solve for
                        # the flux and store the flux in variable x.

			self.x = np.dot(Ainv,source)

			# Calculate the fission source in each spatial bin

			source[:] = self.x[:]*options.delta*2.5*fisXS[:]
                        			
			# Perform the source normalization by dividing by k_eff

			k = sum(source)
                        source[:]=source[:]/k

			# Calculate the difference in the source between
        		# consecutive iterations and take the infinity norm.

			errorDiff[:]=abs(lastSource[:]-source[:])
			error=max(errorDiff)

			# Print statement to show eigenvalue convergence by
                        # iteration.
			#print(j,k,error)

		
                # Calculate the (unnormalized) power of the reactor

                energyPerFission = 200.
		power = self.x[:]*options.delta*energyPerFission*fisXS[:]
                                
		# Print the final output from the diffusion solution.
                
		print "Eigenvalue: ", k, "(", error, ")"
                #print "Max flux: ", max(self.x)[0]
                #print "Total source: ", sum(source)
                #print "Min fisXs: ", min(fisXS)
                #print "Max fixXs: ", max(fisXS)
                #print "Power: ", power[0]
                #print ""

		# Normalize the flux to the appropriate power level.
                # We should probably read the power level from the input file.
                
                powerLevel = 100.
		self.x[:]=self.x[:]*powerLevel/power[0]
		
                #for i in range(0,len(self.x)):
		#	self.x[i]=self.x[i]*100./power[0]
		#sum=0
		#for k in range(1, options.numGroups+1):
		#	for i in range(options.numBins*(k-1),options.numBins*k):
		#		sum=sum+self.x[i]*options.delta/options.numBins*200*NDarray[i-options.numBins*(k-1),0]*Ndata['fuel']['fisxs'][k]
		#print sum
		
		
		
		
