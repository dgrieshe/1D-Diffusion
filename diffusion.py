#last edited by Miriam Rathbun on 6/11/2016
#this script is run from reader.sh. It fills the A matrix which is the linear system of the discrete diffusion equations
#this script requires all the input files to give the same value for "length"



#use sudo apt-get install python.numpy if needed
import numpy as np
import sys, os
import matplotlib.pyplot as plt



#names the variables that came in from the input file
length= float(sys.argv[1])								
#I should make this a condition: if rightBC= 'escape' then apply that diffusion coef rule, eventually
rightBC = sys.argv[2]									
leftBC = sys.argv[3]
numgroups = int(sys.argv[4])
#name of the input file
f = sys.argv[5]											
#retains only the number associated with the input file. Ex: "input1" becomes "1." This was made so the ouput could be nice
name = f[len(f)-5]										 
#number of bins
n = int(sys.argv[6])
#cross sections for each bin. WORK IN PROGRESS
cross_section=np.zeros(n)
#for i in range(0,n):
	#cross_section[i]= float(sys.argv[7+i])
#print cross_section



#this section defines and fills the finite difference coefficient matrix
#note: counting starts at 0 for rows and cols
finite_dif_coef=np.zeros(n+2)	
#cross section . Should eventually give each bin its own cross section. Should be read from input file
xs=[0.75]
#diffusion coefficients. Should be read from input file, eventually. D=1/(3*xs)
diff_coef=np.zeros(n+1)
for i in range(0,n+1):
	diff_coef[i]=1/(3*xs[0])
#width of mesh cell
delta=length/n
#arbitrary source term, constant
sol=[5]


#this section should be incorporated into the section below which fills A
for d in range(0,n+2):
	if d==0: 
		finite_dif_coef[d]=2*diff_coef[d]/(delta+delta*diff_coef[d])
	else:
		if d==n+1:
			finite_dif_coef[d]=2*diff_coef[d-1]/(delta+delta*diff_coef[d-1])
		else:
			finite_dif_coef[d]=2*diff_coef[d-1]*diff_coef[d]/(delta*diff_coef[d-1]+delta*diff_coef[d])


#This section defines and fills matrices A (linear system of diffusion equations) and B (source)
A=np.zeros((n,n)) 										
B=np.zeros((n,1))										
det=np.linalg.det(A) 
#if det=0, A is singular, signal an error, eventually.
while det==0:											
	
   for row in range(0,n):
	   for col in range(0,n):
	   	   if row == col:
	   	   	   if row == 0:
	   	   	   	   A[row,col]=finite_dif_coef[row]+xs[0]*delta
	   	   	   else:
	   	   	   	   if row == n-1:
	   	   	   	   	   A[row,col]=finite_dif_coef[row+2]+xs[0]*delta
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

#calculating the solution x to Ax=B
Ainv=np.linalg.inv(A)
x=np.dot(Ainv,B)
#print A
#print Ainv
#print B
#print x

#plotting solution x
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




