#class created to construct linear system of equations

from material import *

class Construct():
	
###############################################################
	
	def __init__(self):
		self
		
###############################################################
		
	def constructA(self, options, DiffCoef, Gscat, XS, NDarray):
		import numpy as np

		
		nGrps=options.numGroups
		nBins=options.numBins
		d=options.delta
		self.A=np.zeros((nBins*nGrps,nBins*nGrps))


		for k in range(1,nGrps+1):
			for row in range(nBins*(k-1),nBins*k):
				for col in range(nBins*(k-1),nBins*k):
					if row == col:
						if row == nBins*(k-1):
							self.A[row,col]=2*DiffCoef[row]/(d*1+4*DiffCoef[row])+XS[row]*d+2*DiffCoef[row]*DiffCoef[row+1]/(d*DiffCoef[row+1]+d*DiffCoef[row])-Gscat[row-nBins*(k-1),nGrps*(k-1)+k-1]*d
						elif row == nBins*k-1:
							self.A[row,col]=2*DiffCoef[row]/(d*1+4*DiffCoef[row])+XS[row]*d+2*DiffCoef[row-1]*DiffCoef[row]/(d*DiffCoef[row]+d*DiffCoef[row-1])-Gscat[row-nBins*(k-1),nGrps*(k-1)+k-1]*d
						else:
							self.A[row,col]=2*DiffCoef[row]*DiffCoef[row+1]/(d*DiffCoef[row+1]+d*DiffCoef[row])+2*DiffCoef[row-1]*DiffCoef[row]/(d*DiffCoef[row]+d*DiffCoef[row-1])+XS[row]*d-Gscat[row-nBins*(k-1),nGrps*(k-1)+k-1]*d
						#print self.A[row,col]
					elif col == row-1:
						self.A[row,col]=-2*DiffCoef[row-1]*DiffCoef[row]/(d*DiffCoef[row]+d*DiffCoef[row-1])
						#print self.A[row,col]
					elif col == row+1:
						self.A[row,col]=-2*DiffCoef[col-1]*DiffCoef[col]/(d*DiffCoef[col]+d*DiffCoef[col-1])
						#print self.A[row,col]
			for row in range(nBins,nBins*nGrps):
				for col in range(0,nBins*nGrps):
					if row == col+nBins*k:
						i=1
						j=0
						a=row-nBins*k*i
						b=col-nBins*k*j
						while a >= nBins:
							i=i+1
							a=row-nBins*k*i
						while b >= nBins:
							j=j+1
							b=col-nBins*k*j
						self.A[row,col]=-Gscat[a,nGrps*i*k+j]*d

		#print self.A
				
###############################################################
	def invertA(self):
		import numpy as np
		self.inv=np.linalg.inv(self.A)
		
###############################################################
		
		
		
		
		
#Filling matrices A (linear system of diffusion equations) and B (source)

#        # Alternative algorithm for constructing the linear system of equations.
#        # This algorithm offers several potential advantages:
#        #    - Simpler loop indexing scheme
#        #    - Avoids looping over zero-value matrix entries
#        #    - Avoids separate spatial loop for inscattering source
#        #    - Can be used with sparse matrix representation
#        #    - Easy to switch between Marshak and zero flux boundary conditions
#        #    - Source code fits in 120 columns
#        #    - Easy to convert cross section representation, which is important if
#        #      we are computing cross sections on-the-fly.
#
#        nGrps = options.numGroups
#        nBins = options.numBins
#        rank  = nBins*nGrps
#
#        A = np.zeros((rank, rank))
#        B = np.zeros((rank,1))
#
#        for g in range(0, nGrps):
#           for x in range(0, nBins):
#
#              # Determine the unique phase index (cell and energy group) for this entry.
#              # This value is used to index into the cross section arrays.
#
#              i = g*nBins+x
#
#              # Calculate the average diffusion theory between the current cell and the
#              # adjacent cell to the left.  Note that we set coefficients consistent with
#              # the Marshak boundary condition for the cells on the left boundary.  For a
#              # zero flux boundary condition, set deltaAdj equal to zero instead of 4 here.
#
#              if x == 0:
#                 dCoefAdj = 1
#                 deltaAdj = 4     ## Marshak escape boundary condition
#                 #deltaAdj = 0    ## Zero flux boundary condition
#              else:
#                 dCoefAdj = diff_coef[i-1]
#                 deltaAdj = options.delta
#
#              dLeft = (2*diff_coef[i]*dCoefAdj)/(deltaAdj*diff_coef[i] + options.delta*dCoefAdj)
#
#              # Calculate the average diffusion theory between the current cell and the
#              # adjacent cell to the right.  Note that we set coefficients consistent with
#              # the Marshak boundary condition for the cells on the right boundary.  For a
#              # zero flux boundary condition, set deltaAdj equal to zero instead of 4 here.
#
#              if x == (nBins-1):
#                 dCoefAdj = 1
#                 deltaAdj = 4     ## Marshak escape boundary condition
#                 #deltaAdj = 0    ## Zero flux boundary condition
#              else:
#                 dCoefAdj = diff_coef[i+1]
#                 deltaAdj = options.delta
#
#              dRight = (2*diff_coef[i]*dCoefAdj)/(deltaAdj*diff_coef[i] + options.delta*dCoefAdj)
#
#              # Set the matrix coefficient for the cell.  Note that this does not include the
#              # within-group inscattering source term, as this will be included later.
#
#              A[i,i] = dLeft+dRight+xs[i]*options.delta
#              #print(i, i, A[i,i])
#
#              # Set the matrix coefficients for the adjacent cells, except when processing a 
#              # boundary cell.
#
#              if not x == 0:
#                 A[i,i-1] = -dLeft
#                 #print(i, i-1, A[i,i-1])
#
#              if not x == (nBins-1):
#                 A[i,i+1] = -dRight
#                 #print(i, i+1, A[i,i+1])
#
#              # Here we include the group-to-group inscattering source terms.  Because we are
#              # dealing with small 1D systems these terms have been included in the coefficient matrix
#              # rather than lagged in the source term, thus eliminating the need for an iteration
#              # (or sweep) over the inscattering term.
#
#              # We have assumed no upscattering here, but this is not a requirement.  Simply
#              # change the range of the incident neutron groups to include upscattering. 
#
#              for gIn in range(0, g+1):
#                 iIn = gIn*nBins+x
#                 A[i,iIn]=A[i,iIn]-scat[g,gIn]
#                 #print(i, iIn, A[i,iIn])
#
#              # Finally we set the source term for the cell/energy group.
#
#              B[i,0]=source[i]*options.delta
