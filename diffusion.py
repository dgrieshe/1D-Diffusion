#last edited by Miriam Rathbun on 7/18/2016
#This script solves the discrete diffusion equations
#this script can run multiple input files with restrictions



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np    
from nuclide import *
from diffOpts import *
from plotter import * 
from GetFileName import *

fn=FileName()
fn.GetFileName()

for f in fn.listfn:
	print f
	#retains only the number associated with the input file. Ex: "input1" becomes "1" 
	#This was made to make the output name nice
	name = f[len(f)-5]
	
	#variables
	u_g=0.33
	
	options=DiffusionOpts1D()
	options.read(f)
	
	#cross sections for each bin.
	xs=[]
	diff_coef=[]
	source=np.zeros(options.numBins*options.numGroups)
	scat=np.zeros((options.numGroups,options.numGroups))
	
	#creates slab that's all U-235
	for k in range(1,options.numGroups+1):
		M= ('M%i' %k)
		M=Nuclide('U235(%i)' %k)
		M.read()
		for i in range(0,options.numBins):
			xs.append(M.totalXs)
			diff_coef.append(1/(3*(M.totalXs-u_g*M.scatXs)))
		for j in range(0,options.numGroups):
			scat[j,k-1]=eval('M.scatXs%i' %(j+1))
			
			
				
    	#if i<=options.numBins/2:
    	#	M=Nuclide('U235(%i)' % j)
    	#	M.read()
    	#	xs[i]= M.absXs                                  
    	#	diff_coef[i]=1/(3*M.absXs)
    	#else:
    	#	M= Nuclide('U238(%i)' %j)
    	#	M.read()
    	#	xs[i]= M.absXs                                       
    	#	diff_coef[i]=1/(3*M.absXs)       
    		
    #arbitrary source constant source term  		
		
	source[10]=0.1
	#for i in range(0,len(source)):
	#	source[i]=0.1
	
###############################################################################
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

        # Original algorithm for constructing the linear system of equations.

        A=np.zeros((options.numBins*options.numGroups,options.numBins*options.numGroups))
        B=np.zeros((options.numBins*options.numGroups,1))

        for k in range(1,options.numGroups+1):
                for row in range(options.numBins*(k-1),options.numBins*k):
                        for col in range(options.numBins*(k-1),options.numBins*k):
                                if row == col:
                                        if row == options.numBins*(k-1):
                                                A[row,col]=2*diff_coef[row]/(options.delta*1+4*diff_coef[row])+xs[row]*options.delta+2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])-scat[k-1,k-1]
                                                #print(row, col, A[row,col])
                                        elif row == options.numBins*k-1:
                                                A[row,col]=2*diff_coef[row]/(options.delta*1+4*diff_coef[row])+xs[row]*options.delta+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])-scat[k-1,k-1]
                                                #print(row, col, A[row,col])
                                        else:
                                                A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(options.delta*diff_coef[row+1]+options.delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])+xs[row]*options.delta-scat[k-1,k-1]
                                                #print(row, col, A[row,col])
                                elif col == row-1:
                                        A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(options.delta*diff_coef[row]+options.delta*diff_coef[row-1])
                                        #print(row, col, A[row,col])
                                elif col == row+1:
                                        A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(options.delta*diff_coef[col]+options.delta*diff_coef[col-1]) 
                                        #print(row, col, A[row,col])
                for row in range(options.numBins,options.numBins*options.numGroups):
                        for col in range(0,options.numBins*options.numGroups):
                                if row == col+options.numBins*k:
                                        A[row,col]=-scat[k,k-1]
                                        #print(row, col, A[row,col])
        for row in range(0,len(source)):
                B[row,0]=source[row]*options.delta

###############################################################################   

    
     
###############################################################################
#calculating the solution x to Ax=B   
	Ainv=np.linalg.inv(A)
	x=np.dot(Ainv,B)
	print x[0]

	
	results=Plotter()
	results.plot(x,1,options.numBins,options.numGroups,name)
    
###############################################################################
