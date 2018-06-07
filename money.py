

# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys


#lambda 2, phi contour boundary, gamma perscription on top with Vknob
rhoEdge1 = 0.5625
rhoC1  = array([1.0, 1.2054, 1.332, 1.49, 1.6913, 1.931, 2.145, 2.438, 2.766, 3.121, 3.6131, 4.2527, 4.8065 ])
mEx1   = array([0.0, 0.041, 0.082, 0.123, 0.164, 0.205, 0.246, 0.287, 0.328, 0.369, 0.410, 0.451, 0.492 ])
Vknob1 = array([0.0, -0.05, -0.1, -0.15, -0.15, -0.15, -0.2, -0.2, -0.2, -0.21, -0.2, -0.2, -0.2 ])
rhoB1 = rhoC1/rhoEdge1
fEx1 = (1.0/0.82348)*mEx1
print(rhoB1)

#lambda 2, phi contour boundary, gamma perscription on top with Vknob
rhoEdge = 0.25
rhoC  = array([1.0, 1.2030, 1.3656, 1.55010, 1.8201, 1.973, 2.161, 2.385, 2.826,  3.0872, 3.286, 3.491])
mEx   = array([0.0, 0.05,   0.10,   0.15,    0.20,   0.225, 0.25,  0.275, 0.28125, 0.2875, 0.29375, 0.30])
Vknob = array([0.0, 0.08,   0.13,   0.177,   0.26,   0.30,  0.35,  0.40,  0.50,   0.55,   0.58,  0.61])
rhoB = rhoC/rhoEdge
fEx = (1.0/1.6469)*mEx


plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)
plt.axis([1.0, 20.0, 0, 0.7])

plt.ylabel(r'$f_{\rm ex}$',fontsize=15)
plt.xlabel(r'$\rho_c/\rho_b$',fontsize=15)

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=10)
ax.tick_params(which='minor',size=7)

plt.plot(rhoB,fEx,c='black',alpha=0.75,label=r'$\lambda=2$')
plt.scatter(rhoB,fEx,c=Vknob)
plt.plot(rhoB1,fEx1,c='blue',alpha=0.75,label=r'$\lambda=1$')
plt.scatter(rhoB1,fEx1,c=Vknob1)
plt.colorbar(label=r'Vknob | $\gamma$')
plt.scatter([3.05459/rhoEdge],[0.28125/1.6469],c='black',s=100,alpha=0.15)

plt.legend(fontsize=10)

'''
ax2 = ax.twiny()
ax2.set_xlabel(r'$\rho_c$',fontsize=15)
ax2.set_xlim( ax.get_xlim() )
ax2.set_xticks( ax.get_xticks() )
ax2.set_xticklabels( around( (rhoEdge)*ax.get_xticks(), decimals=1) )
ax2.tick_params(which='both',width=2)
ax2.tick_params(which='major',size=10)
ax2.tick_params(which='minor',size=7)
'''


plt.savefig('MoneyPlot.png')
