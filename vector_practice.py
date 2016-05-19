import numpy as np
import random

#define a square matrix A
#note: counting starts at 0 for rows and cols
n=random.randint(2,10)
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
print('Solution to Ax=B where A=') 
print(A)
print('and B=')
print(B)
print('is:')
print(x)

#print(B)
#print(np.dot(A,x))										#check that A*x=B


	   




