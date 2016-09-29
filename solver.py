#class created to solve Ax=B

class Solve:
	
	def __init__(self):
		self
		
	def solve(self, options, Ainv, source):
		
		import numpy as np
		B=np.zeros((options.numBins*options.numGroups,1))
		lastSource=np.zeros(options.numBins*options.numGroups)
		error=1000 
		j=0
		while error > 10 ** (-4):
			k=0
			errorDiff=np.zeros(options.numBins*options.numGroups)
			j=j+1
			for i in range(0,len(source)):
				lastSource[i]=source[i]
				B[i]=source[i]*options.delta
			self.x=np.dot(Ainv,B)
			for i in range(0,len(source)):
				source[i]=self.x[i]*2.5*0.006
				k=k+source[i]*options.delta
			for i in range(0,len(source)):
				source[i]=source[i]/k
				errorDiff[i]=abs(lastSource[i]-source[i])
			error=max(errorDiff)
			print(j,k,error)
		
		#One way to avoid looping through each value in the array
		#doesn't seem any faster, because this does seem to be a loop
		#result = [v * 2.5*0.006 for v in x]
		
		
		