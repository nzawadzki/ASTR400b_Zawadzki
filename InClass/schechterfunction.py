#InClass worksheet 1
#1/25/18
#plotting the schechter function

import numpy as np
import matplotlib.pyplot as plt

#constants
h = 70.4/100
Mstar = -23.19-5*np.log(h)
alpha = -0.81
phi_star = 1.66e-2

def Schechter(alpha, Mstar, phi_star, M):
    #define helpful variables to break up function
    a = 10**(0.4*(Mstar-M)*(alpha+1)) #power law
    exp = np.exp(-10**(0.4*(Mstar-M))) #decaying exponential
    Phi = (0.4*np.log(10))*phi_star*a*exp
    return Phi

#etting axis limits
M=np.arange(-26,-17,0.1)
func = Schechter(alpha, Mstar, phi_star, M)
#set numbers to go high to low
plt.xlim(-17,-26)
#set axis to be log based
plt.semilogy(M,func)
plt.show()



