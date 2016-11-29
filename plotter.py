#class to plot


class Plotter:
	
	
	def __init__(self):
		self
		
	def plot(self, solution, group, numBins, numGroups, name, n):
		
		import numpy as np
		import matplotlib.pyplot as plt
		plt.clf()
		groups=[]
		for k in range(1,numGroups+1):
			l=0
			x_group= np.zeros(numBins)
			for i in range(numBins*(k-1),numBins*k):
				x_group[l]=solution[i]
				l=l+1
				
			n=str(n)
			groups.append('group%i' %k)
			colors=['r','b','g','o','n','p']	
			groups[k-1], = plt.plot(x_group, colors[k-1], label=groups[k-1])
			#plt.ylim([0,0.5])
			plt.ylabel('flux')
			plt.legend(groups)
			plt.savefig('./output/figure'+name+n)
	
	