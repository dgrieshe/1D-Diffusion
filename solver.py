#class created to solve Ax=B

class Solve:
	
	def __init__(self):
		self
		
	def solve(self, options, Ainv, source, NDarray, Ndata):
		
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
		norm=0
		for k in range(1, options.numGroups+1):
			for i in range(options.numBins*(k-1),options.numBins*k):
				norm=norm+self.x[i]*options.delta/options.numBins*200*NDarray[i-options.numBins*(k-1),0]*Ndata['fuel']['fisxs'][k]
		for i in range(0,len(self.x)):
			self.x[i]=self.x[i]/norm
		#sum=0
		#for k in range(1, options.numGroups+1):
		#	for i in range(options.numBins*(k-1),options.numBins*k):
		#		sum=sum+self.x[i]*options.delta/options.numBins*200*NDarray[i-options.numBins*(k-1),0]*Ndata['fuel']['fisxs'][k]
		#print sum
		
		
		
		