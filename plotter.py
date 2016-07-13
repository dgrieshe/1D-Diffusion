#class to plot


class Plotter:
	
	
	def __init__(self):
		self
		
	def plot(self, solution, group, numBins, numGroups, name):
		
		import numpy as np
		import matplotlib.pyplot as plt
		plt.clf()
		for k in range(1,numGroups+1):
			l=0
			x_group= np.zeros(numBins)
			for i in range(numBins*(k-1),numBins*k):
				x_group[l]=solution[i]
				l=l+1
				
			plt.plot(x_group)
			plt.ylabel('flux')
			plt.savefig('./output/figure'+name)
	
	