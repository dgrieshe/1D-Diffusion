#last edited by Miriam Rathbun on 5/25/2016
#this script is run from reader.sh. It creates two matrices A and B and solves for x such that Ax=B. Output is dumped to files in the output folder. 

#use sudo apt-get install python.numpy if needed for any of these (required the first time on a new machine)
import numpy as np
import random
import sys, os
import test
import matplotlib.pyplot as plt

#names the variables that came in from the input file
length= int(sys.argv[1])
rightBC = sys.argv[2]
leftBC = sys.argv[3]
numgroups = int(sys.argv[4])
f = sys.argv[5]											#name of the input file
name = f[len(f)-5]										#retains only the number associated with the input file. Ex: "input1" becomes "1." This was made so the ouput could have a pretty name 

#define a square matrix A
#note: counting starts at 0 for rows and cols
n=length  												#just to show that the passed variables can be used in this script
A=np.zeros((n,n)) 										#creates a random square matrix A (dimensions (n,n)) for the Ax=B calc
B=np.zeros((n,1))										#creates vector B of length n

det=np.linalg.det(A)  									#determinant of A
while det==0:											#if det=0, A is singular and should be re-made
	
   for row in range(0,n):								#fill A and B with random numbers between 0 and 10
	   for col in range(0,n):
		   A[row,col]=random.uniform(0,10)
	   B[row,0]=random.uniform(0,10)				    
   det=np.linalg.det(A)


Ainv=np.linalg.inv(A)									#take the inverse of A
x=np.dot(Ainv,B)										#solution x= matrix multiplication of Ainv and B


f = open('./output/output_'+name+'.text', 'w')
print >> f, 'Solution to Ax=B where A= \n', A
print >> f, 'and B= \n', B
print >> f, 'is: \n', x
f.close()

#print(B)
#print(np.dot(A,x))										#check that A*x=B




