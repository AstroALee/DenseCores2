# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

cmapalpha = 0.0

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

plt.axis([0, rPlotmax, 0, zPlotmax])

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=5 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

#moves numbers out a little
ax.tick_params(axis='y',pad=9)
ax.tick_params(axis='x',pad=9)

# Max and Min for ColorMap
# Overwrite min max plotted
#plotmin = 0
#plotmax = 20
cmapmin = plotmin
cmapmax = plotmax

# ColorMap choice
cmap = 'Paired'
cmap = 'RdBu_r'
#cmap = 'Spectral'
#cmap = 'Blues_r'
if(cmapalpha>0):
    plt.pcolor(RGrid,ZGrid,PreparedPlot,cmap=cmap,vmin=cmapmin,vmax=cmapmax,alpha=cmapalpha)
    plt.colorbar()


# override
CLevels = [ 0.5, 0.75, 1, 1.25, 1.5  ]
CS = plt.contour(RMesh,ZMesh,PreparedPlot,levels=CLevels,colors='k',linewidths=2.5)
#plt.clabel(CS, inline=1, fontsize=9)

# Filament Boundary
area = pi*2.0
#plt.scatter(VCont,(ZGrid[:-1]+0.5*DeltaZ),s=1.5*area,c='black')
plt.plot(VCont,(ZGrid[:-1]+0.5*DeltaZ),c='black',linewidth=2.5,alpha=0.3)

# B fields
if(Bfield=='Y'):
    for i in range(NBCont):
        plt.plot(BContours[i],(ZGrid[:-1]+0.5*DeltaZ),'--',c='black',linewidth=2,alpha=0.7)

#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)


# Creates Plot
#plt.show()
plt.savefig('Contour.pdf')
