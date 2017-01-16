#class to plot


class Plotter:
	
####################################################################
	def __init__(self):
		self

####################################################################		
	def plotFLUX(self, solution, group, numBins, numGroups, inp, n):
		
		import numpy as np
		import matplotlib.pyplot as plt
		plt.clf()
		groups = []
		for k in range(1,numGroups+1):
			l = 0
			x_group = np.zeros(numBins)
			for i in range(numBins*(k-1),numBins*k):
				x_group[l] = solution[i]
				l = l+1
				
			n = str(n)
			inp = str(inp)
			groups.append('group%i' %k)
			colors = ['r','b','g','o','n','p']	
			groups[k-1], = plt.plot(x_group, colors[k-1], label=groups[k-1])
			plt.ylabel('flux')
			plt.legend(groups)
			plt.savefig('./output/flux'+inp+n)


####################################################################
	def plotFluxRMS(self, sol, inp, n, flux1, flux2):

		 import numpy as np
		 import matplotlib.pyplot as plt
		 plt.clf()

		 RMS = []
		 for i in range(0,len(flux1)):
		 	RMS.append(np.sqrt(sum((flux1[i][:]-flux2[i][:])**2)))

		 plt.plot(RMS)
		 plt.savefig('./output/FluxRMS')





####################################################################
	def plotMASSU(self, massU, n, inp):

		import numpy as np
		import matplotlib.pyplot as plt
		plt.clf()
		mass = np.zeros(n)
		legend = np.zeros(inp)
		colors = ['r','b','g','o','n','p']	

		for i in range(0,inp):
			for j in range(0,n):
				mass[j] = massU[n*i+j]
				#if i == 0 and j == 0:
				#	print mass[j]
			legend[i] = i+1
			plt.plot(mass, colors[i])
			plt.legend(legend)
		#print mass[j]
		plt.savefig('./output/massU')


		#n = str(n)
		#plt.clf()
		#plt.plot(massU)
		#plt.savefig('./output/massU'+name+n)


		



	
	