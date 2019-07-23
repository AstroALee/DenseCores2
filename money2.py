

# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys


#lambda 1, phi contour boundary, gamma perscription on top with Vknob (4.8065 last real entry)
rhoEdge1 = 0.5625
rhoC1  = array([1.0, 1.2054, 1.332, 1.49, 1.6913, 1.931, 2.145, 2.438, 2.766, 3.121, 3.6131, 4.2527, 4.8065, 5.3, 5.9, 6.5,7,8,9,10,10.5 ])
mEx1   = array([0.0, 0.041, 0.082, 0.123, 0.164, 0.205, 0.246, 0.287, 0.328, 0.369, 0.410, 0.451, 0.492,0.533, 0.574, 0.615,0.656, 0.687,0.7,0.705,0.71 ])
Vknob1 = array([0.0, -0.05, -0.1, -0.15, -0.15, -0.15, -0.2, -0.2, -0.2, -0.21, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2 , -0.2, -0.2,-0.25,-0.25,-0.25 ])
rhoB1 = rhoC1/rhoEdge1
fEx1 = (1.0/0.82348)*mEx1
print(rhoB1)

#lambda 2, phi contour boundary, gamma perscription on top with Vknob
rhoEdge = 0.25
rhoC  = array([1.0, 1.2030, 1.3656, 1.55010, 1.8201, 1.973, 2.161, 2.385, 2.826,  3.0872, 3.286, 3.491])
mEx   = array([0.0, 0.05,   0.10,   0.15,    0.20,   0.225, 0.25,  0.26, 0.265, 0.28125, 0.2875, 0.29375])
Vknob = array([0.0, 0.08,   0.13,   0.177,   0.26,   0.30,  0.35,  0.40,  0.50,   0.55,   0.58,  0.61])
rhoB = rhoC/rhoEdge
fEx = (1.0/1.6469)*mEx

# beta 0.1, lambda 1
rhoEdge12 = 0.5625
rhoC12  = array([1.0,   1.49, 1.6913,   2.821, 3.6131, 4.2527,  5.3,  6.5,   9 ,10.5, 13,15, 17])
mEx12   = array([0.0,   0.123, 0.164,   0.37,0.45, 0.52,0.59 , 0.68,0.8,0.84, 0.86,0.83,0.8])
Vknob12 = array([0.0, -0.05, -0.1, -0.15, -0.15, -0.15, -0.2 , -0.2, -0.2, -0.2])
rhoB12 = rhoC12/rhoEdge12
fEx12 = (1.0/0.82348)*mEx12


plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)
plt.axis([1.0, 40.0, 0, 1.4])

plt.ylabel(r'$f_{\rm ex}$',fontsize=15)
plt.xlabel(r'$\rho_c/\rho_b$',fontsize=15)

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=10)
ax.tick_params(which='minor',size=7)

#plt.plot(rhoB,fEx,c='black',alpha=0.75,label=r'$\lambda=2$')
#plt.scatter(rhoB,fEx,c=Vknob)
#plt.plot(rhoB1,fEx1,c='blue',alpha=0.35,label=r'$\lambda=1$')
#plt.scatter(rhoB1,fEx1,c=Vknob1)


#plt.scatter([12.8/rhoEdge1],[0.82],c='black',s=50,alpha=0.45)
size=80
#plt.scatter([12.8/rhoEdge1],[0.81],c='blue',s=size,alpha=0.45)
#plt.scatter([13.65/rhoEdge1],[0.77],c='blue',s=size,alpha=0.45)
#plt.scatter([14.65/rhoEdge1],[0.72],c='blue',s=size,alpha=0.45)
#plt.scatter([14.65/rhoEdge1],[0.7],c='black',s=50,alpha=0.15)
#plt.scatter([15.15/rhoEdge1],[0.65],c='black',s=50,alpha=0.45)


#rhoC1_2 = array([1.0,     1.49, 1.6913,  2.438, 2.766,  3.6131, 4.2527,  5.3,  6.5, 7.5,   9 ,10.5, 12])
#mEx1_2   = array([0.0,   0.123, 0.164,   0.287, 0.328,  0.410, 0.459,0.538 , 0.625, 0.665, 0.69, 0.7, 0.69 ])
#rhoC1_2 = array([1.0,  1.2,   1.49, 1.6913, 2.2, 2.438, 2.766,  3.2131])
#mEx1_2   = array([0.0, 0.07,  0.123, 0.164, 0.24,  0.287, 0.328,  0.40  ])


rhoC1_2 = array([1.0,     1.49, 1.6913,  2.438, 2.766,  3.6131, 4.2527,  5.3,  6.5, 7.5,   9 ,10.5, 12,13.5,15])
mEx1_2   = array([0.0,   0.123, 0.164,   0.287, 0.328,  0.410, 0.459,0.538 , 0.625, 0.665, 0.69, 0.7, 0.69,0.67,0.625 ])
plt.scatter(0.9*rhoC1_2/rhoEdge1,(1.2/0.82348)*mEx1_2,c='blue',s=80,alpha=0.45,label=r'$\lambda=1$ high-res, $\beta=1$')


rhoC12  = array([1.0,   1.49, 1.6913, 2.25,  2.821, 3.6131, 4.2527,  5.3,  6.5, 7.8,  9 ,10.5, 13,15, 17])
mEx12   = array([0.0,   0.123, 0.164, 0.25,  0.37,0.45, 0.52,0.59 , 0.68,       0.74, 0.8,0.84, 0.86,0.85,0.82])
plt.scatter(0.8*rhoC12/rhoEdge12,(1.2/0.82348)*mEx12,c='black',s=80,alpha=0.45,label=r'$\lambda=1$, high-res, $\beta=0.1$')



#plt.scatter(rhoC1_2/rhoEdge1,(1.0/0.82348)*mEx1_2,c='blue',s=80,alpha=0.45,label=r'$\lambda=1$ high-res, $\beta=1$')
#plt.scatter(rhoC1_2/rhoEdge1,1.5*(1.0/0.82348)*mEx1_2,c='blue',s=80,alpha=0.45,label=r'$\lambda=1$ high-res, $\beta=1$')

#plt.scatter(rhoC12/rhoEdge12,(1.0/0.82348)*mEx12,c='black',s=80,alpha=0.45,label=r'$\lambda=1$, high-res, $\beta=0.1$')
plt.title(r'$n=2/3$')

#plt.colorbar(Vknob,label=r'Vknob | $\gamma$')
#plt.scatter([3.05459/rhoEdge],[0.28125/1.6469],c='black',s=70,alpha=0.15)

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


plt.savefig('MoneyPlot2.png')
