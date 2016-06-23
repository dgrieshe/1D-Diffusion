#last edited by Miriam Rathbun on 6/22/2016
#this script is run from reader.sh. It fills the A matrix which is the linear system of the discrete diffusion equations
#this script can run multiple input files with restrictions



#imports
#use sudo apt-get install python.numpy if needed
import numpy as np
import sys, os
import matplotlib.pyplot as plt



#names the variables that came in from the input file
length= float(sys.argv[1])								
#Eventually, I should make this a condition: if rightBC= 'escape' then apply that diffusion coef rule
rightBC = sys.argv[2]									
leftBC = sys.argv[3]
numgroups = int(sys.argv[4])
#name of the input file
f = sys.argv[5]											
#retains only the number associated with the input file. Ex: "input1" becomes "1." This was made to make output name nice
name = f[len(f)-5]										 
#number of bins
n = int(sys.argv[6])
#cross sections for each bin. should denote an error when number of xs is different from n
xs=np.zeros(n)
for i in range(0,n):
	xs[i]= float(sys.argv[7+i])
#print cross_section
#diffusion coefficients. Should be read from input file, eventually. D=1/(3*xs)
diff_coef=np.zeros(n+1)
for i in range(0,n):
	diff_coef[i]=1/(3*xs[i])
#width of mesh cell
delta=length/n
#arbitrary source term, constant
source=[5]



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
	   	   	   	   #this uses the ghost bin. I don't know what to assign it yet, so I'm giving it a value of 1. 
	   	   	   	   A[row,col]=2*diff_coef[row]*1/(delta*1/+delta*diff_coef[row])+xs[row]*delta
	   	   	   else:
	   	   	   	   if row == n-1:
	   	   	   	   	   #this uses the ghost bin. I don't know what to assign it yet, so I'm giving it a value of 1.
	   	   	   	   	   A[row,col]=2*diff_coef[row]*1/(delta*1/+delta*diff_coef[row])+xs[row]*delta
	   	   	   	   else:
	   	   	   	       A[row,col]=2*diff_coef[row]*diff_coef[row+1]/(delta*diff_coef[row+1]+delta*diff_coef[row])+2*diff_coef[row-1]*diff_coef[row]/(delta*diff_coef[row]+delta*diff_coef[row-1])+xs[row]*delta
	   	   else:
	   	   	   if col == row-1:
	   	   	   	   A[row,col]=-2*diff_coef[row-1]*diff_coef[row]/(delta*diff_coef[row]+delta*diff_coef[row-1])
	   	   	   else:
	   	   	   	   if col == row+1:
	   	   	   	   	   A[row,col]=-2*diff_coef[col-1]*diff_coef[col]/(delta*diff_coef[col]+delta*diff_coef[col-1]) 
	   B[row,0]=source[0]*delta		    
   det=np.linalg.det(A)
   
   
#calculating the solution x to Ax=B
Ainv=np.linalg.inv(A)
x=np.dot(Ainv,B)



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



#This section experiments with dictionaries
#I want to make the dictionary have a better/more complex structure
#I also want to read the dictionary items in from input
nuclides = {}
nuclides['U-235']=[2, 3, 5]
nuclides['U-238']=[4, 5, 5]                   
nuclides['U-239']=[6, 7, 1]   

nuclide_choice = raw_input("What nuclide would you like to know about? Choices are: U-235, U-238, U-239 : ")

def read_nuclide_type(nuclides, nuclide_choice):
	for nuclide, nuclide_type in nuclides.items():
		if nuclide == nuclide_choice:
			print("The nuclide %s has:\n"
				"scattering xs = %s \n"
				"absorbtion xs = %s \n"
				"fission xs = %s" % (nuclide, nuclide_type[0], nuclide_type[1], nuclide_type[2]))

read_nuclide_type(nuclides, nuclide_choice)

#del nuclides['U-239']	
#print("The nuclide %s " % 'U-235')
#print("is %s" % nuclides['U-235'][0])

