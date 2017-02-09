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

        plt.clf()


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

         plt.clf()


####################################################################
    def plotINTflux(self, time, powerPlot, inp):

        import numpy as np
        import matplotlib.pyplot as plt
        import csv

        #print powerPlot

        #legend = np.arange(1,5)
        #colors = ['r','b','g', 'c', 'm', 'y', 'k']  

        #plt.plot(powerPlot, colors[inp])
        #plt.legend(legend)
        #plt.savefig('./output/INTflux')


        if inp == 1:
            f = open('output/filename', 'w')
            f.close()

        f = open('output/filename', 'a')
        for i in range(0,len(powerPlot)):
            f.write(str(time[i]))
            f.write(', ')
        f.write('\n')
        for i in range(0,len(powerPlot)):
            f.write(str(powerPlot[i]))
            f.write(', ')
        f.write('\n')
        f.write('\n')
        # I need it to append to the end of the file each time.




####################################################################
    def plotMASSU(self, timeMU, massU, massU1, n):

        import numpy as np
        import matplotlib.pyplot as plt
        plt.clf()
        mass = np.zeros(len(timeMU))
        #legend = np.zeros(inp)
        #colors = ['r','b','g', 'c', 'm', 'y', 'k']

        for j in range(0,len(timeMU)):
            mass[j] = massU1-massU[j]
                #if i == 0 and j == 0:
                #   print mass[j]

        f = open('output/massU', 'a')
        for k in range(0,len(timeMU)):
            f.write(str(timeMU[k]))
            f.write(', ')
        f.write('\n')
        for k in range(0,len(timeMU)):
            f.write(str(mass[k]))
            f.write(', ')
        f.write('\n')
        f.write('\n')


        #legend[i] = i+1
        #plt.plot(mass, colors[i])
        #plt.legend(legend)
        #print mass[j]
        #plt.savefig('./output/massU')


