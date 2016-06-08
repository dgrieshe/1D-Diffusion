#last edited by Miriam Rathbun on 5/26/2016
#this script is run from reader.sh. It creates two matrices A and B and solves for x such that Ax=B. Output is dumped to files in the output folder. 
#this script requires all the input files to give the same value for "length"

#use sudo apt-get install python.numpy if needed for any of these (required the first time on a new machine)
import numpy as np
import sys, os
import matplotlib.pyplot as plt

#names the variables that came in from the input file
length= float(sys.argv[1])								#length of the slab. For input1= 400cm
rightBC = sys.argv[2]									#I should make this mean that if rightBC= 'escape' then apply that diffusion coef 
leftBC = sys.argv[3]
numgroups = int(sys.argv[4])
f = sys.argv[5]											#name of the input file
name = f[len(f)-5]										#retains only the number associated with the input file. Ex: "input1" becomes "1." This was made so the ouput could have a pretty name 
n = int(sys.argv[6])									#number of bins

#define a square matrices
#note: counting starts at 0 for rows and cols
finite_dif_coef=np.zeros(n+2)							#defines the matrix that will be populated by the finite difference coupling coefficients, d
diff_coef=[1,2,3,4,5]									#data of diffusion coefficients for each boundary. To be filled with real values. This should be read from the input file
xs=[3000] 												#cross section data
delta=length/n											#defines the width of the mesh cell
sol=[5]													#arbitrary source term, constant
for d in range(0,n+2):									#populates finite_dif_coef using equations derived from Diffusion-Notes_Rev1 and data from diff_coef
	if d==0: 
		finite_dif_coef[d]=2*diff_coef[d]/(delta+delta*diff_coef[d])				#depends on the fact that deltas are constant and assumes that D^-1=1
	else:
		if d==n+1:
			finite_dif_coef[d]=2*diff_coef[d-1]/(delta+delta*diff_coef[d-1])		#depends on the fact that deltas are constant and assumes that D^5=1
		else:
			finite_dif_coef[d]=2*diff_coef[d-1]*diff_coef[d]/(delta*diff_coef[d-1]+delta*diff_coef[d])	#depends on deltas constant

A=np.zeros((n,n)) 										#creates a random square matrix A (dimensions (n,n)) for the Ax=B calc
B=np.zeros((n,1))										#creates vector B of length n


det=np.linalg.det(A)  									#determinant of A
while det==0:											#if det=0, A is singular and should be re-made
	
   for row in range(0,n):								#fill A and B with random numbers between 0 and 10
	   for col in range(0,n):
	   	   if row == col:
	   	   	   if row == 0 or row == n-1:
	   	   	   	   A[row,col]=finite_dif_coef[row]+xs[0]*delta
	   	   	   else:
	   	   	   	   A[row,col]=finite_dif_coef[row+2]+finite_dif_coef[row+1]+xs[0]*delta
	   	   else:
	   	   	   if col == row-1:
	   	   	   	   A[row,col]=-finite_dif_coef[col+2]
	   	   	   else:
	   	   	   	   if col == row+1:
	   	   	   	   	   A[row,col]=-finite_dif_coef[row+2]
	   B[row,0]=sol[0]*delta		    
   det=np.linalg.det(A)


Ainv=np.linalg.inv(A)									#take the inverse of A
x=np.dot(Ainv,B)										#solution x= matrix multiplication of Ainv and B

plt.plot(x)
plt.ylabel('some numbers')
plt.savefig('./output/figure'+name)

if name=='1':
    f = open('./output/output.text', 'w')
    f.close()

f = open('./output/output.text', 'a')

i=0
while i < n:
   print >> f, float(x[i])
   i = i+1
f.close()

#print(B)
#print(np.dot(A,x))										#check that A*x=B




