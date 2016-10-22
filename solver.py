#class created to solve Ax=B

class Solve:
	
	def __init__(self):
		self
		
	def solve(self, options, Ainv, source, fisXS, NDarray, Ndata):
		
		import numpy as np
		B=np.zeros((options.numBins*options.numGroups,1))
		lastSource=np.zeros(options.numBins*options.numGroups)
		error=1000 
		j=0
		while error > 10 ** (-5):
			k=0
			errorDiff=np.zeros(options.numBins*options.numGroups)
			j=j+1
			for i in range(0,len(source)):
				lastSource[i]=source[i]
				B[i]=source[i]
			self.x=np.dot(Ainv,B)

			for i in range(0,len(source)):
				source[i]=self.x[i]*options.delta*2.5*fisXS[i]
				k=k+source[i]

			for i in range(0,len(source)):
				source[i]=source[i]/k
				errorDiff[i]=abs(lastSource[i]-source[i])
			error=max(errorDiff)
			#print(j,k,error)

		power=0
                for i in range(0,len(self.x)):
			power=self.x[i]*options.delta*200*fisXS[i]
                
		print "Eigenvalue: ", k, "(", error, ")"
                #print "Max flux: ", max(self.x)[0]
                #print "Total source: ", sum(source)
                #print "Min fisXs: ", min(fisXS)
                #print "Max fixXs: ", max(fisXS)
                #print "Power: ", power[0]
                #print ""
		
                for i in range(0,len(self.x)):
			self.x[i]=self.x[i]*100./power[0]
		#sum=0
		#for k in range(1, options.numGroups+1):
		#	for i in range(options.numBins*(k-1),options.numBins*k):
		#		sum=sum+self.x[i]*options.delta/options.numBins*200*NDarray[i-options.numBins*(k-1),0]*Ndata['fuel']['fisxs'][k]
		#print sum
		
		
		
		
