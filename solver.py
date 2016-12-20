# class created to solve Ax=B

class Solve:
	
	def __init__(self):
		self
		
	def solve(self, options, Ainv, source, fisXS, NDarray, Ndata, n):
		
		import numpy as np

		# Create local array variables.
		lastSource = np.zeros(options.numBins*options.numGroups)
		errorDiff = np.zeros(options.numBins*options.numGroups)

		# Initialize local variables
		error = 1000 
		j = 0
                
                # Loop until the source difference between iterations falls
                # below the specified threshold.
                
		while error > options.ConvError:
                	# Reset local variables
			k = 0
			j = j+1

			errorDiff[:] = 0.

			# Save the previous source so that we can compute the
			# residual error for the iteration.
			lastSource[:] = source[:]
                        
			# Take the dot product of A-inverse and B to solve for
            # the flux and store the flux in variable x.

			self.x = np.dot(Ainv,source[:]*options.delta)
			#print sum(source[:]*options.delta)

			# Calculate the fission source in each spatial bin

			source[:] = self.x[:]*options.nYield*fisXS[:]
                        			
			# Perform the source normalization by dividing by k_eff

			k = sum(source)*options.delta
                        source[:] = source[:]/k

			# Calculate the relative difference in the source between
        	# consecutive iterations and take the infinity norm.

			errorDiff[:] = abs((lastSource[:]-source[:])/lastSource[:])
			error = max(errorDiff)

			# Print statement to show eigenvalue convergence by
            # iteration.
			#print(j,k,error)

		#print self.x
		#print source


                                
		# Print the final output from the diffusion solution.
                
		print n, "Eigenvalue: ", k, "(", error, ")"
                #print "Max flux: ", max(self.x)[0]
                #print "Total source: ", sum(source)
                #print "Min fisXs: ", min(fisXS)
                #print "Max fixXs: ", max(fisXS)
                #print "Unnormalized power: ", power[0]
                #print ""		
		
		
		
