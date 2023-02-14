# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:35:04 2023

@author: ellio
"""

import numpy as np
import matplotlib.pyplot as plt

#-----------------------READ ME------------------------#

#To use script you will only need to change SYSTEM VARIABLES for a nitrogen system
#Iterate to find the desired results
#INTIAL PRESSURE will be dictated by what the system can handle and what you can get a hold of at EuRoC
#Iterate VOLUME to achieve the desired total impulse, total impulse is equal to momentum (nozzle radius does not effect total impulse)
#Adjust NOZZLE RADIUS to find achieve different thrust curves, smaller nozzle leads to less thrust but a longer burn and vice versa
#If the code breaks adjust the time range, you want to just see the asymtote of the thrust curve form

#-------------------------------------------------------#
#------------------System Variables---------------------#
#-------------------------------------------------------#


Volume = 0.15            #volume of the pressure vessel in [L]
IntialPressure = 120     #Starting pressure of the pressure vessel in [Bar]
NozzleRadius = 0.5       #Radius of the nozzles throat in [mm]
t = np.arange(0,6,0.01)  #if the script breaks reduce the total time, if you cant see the asymptote increase total time


#-------------------------------------------------------#
#------------------System Constants---------------------#
#-------------------------------------------------------#


gamma = 1.4           #specific heat ratio of the gas being used (1.4 for nitrogen)
R = 300               #gas constant?
T_0_degC = 20         #Initial temperature of the pressure vessel [C]
Cd = 0.6              #Nozzle performance parameter; poor conical nozzle (0.6), perfect bell curve (0.98), frictionless nozzle (1)


#-------------------------------------------------------#
#----------------------Main Script----------------------#
#-------------------------------------------------------#


#convert to SI units
V = Volume/1000
P0 = IntialPressure*10**5
T0 = T_0_degC + 273.15
rho0 = P0/(R*T0)

#Nozzle area, radius in mm 
A = np.pi * (NozzleRadius*10**-3)**2

#Calculate initial speed of sound in the gas
c0 = np.sqrt(gamma*R*T0)

#Calculate tau (constant coefficient defining tank parameters)
tau = (V/(Cd*A*c0))*((gamma+1)/2)**((gamma+1)/(2*(gamma-1)))

#Calculate pressure
P = P0*(1+((gamma-1)/2)*(t/tau))**(2*gamma/(1-gamma))

#Calculate density
rho = rho0*(1+((gamma-1)/2)*(t/tau))**(2/(1-gamma))

#calculate tempertaure
T = T0*(1+((gamma-1)/2)*(t/tau))**-2

#calculate mass flow rate
m_dot = Cd*A*rho*np.sqrt(gamma*R*T)*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))

#calculate thrust
F  = m_dot * np.sqrt(((10**5/P)**((gamma-1)/(-gamma))-1)/((gamma-1)/2)) * c0

#calculate isp
isp = F / (m_dot*9.81)

#nummerical integration of the thrust curve to find total impulse
def trapz(x,y):
    I = 0
    for i in range((len(x))-1):
        I = I + ((y[i]+y[i+1])/2)*(x[i+1]-x[i])
    return(I)

totalImpulse = trapz(t, F)


#-------------------------------------------------------#
#---------------------Display Data----------------------#
#-------------------------------------------------------#


print('Total Impulse [Ns]:',round(totalImpulse,2))
print('Peak Thrust [N]:',round(F[0],2))
print('Peak Specific Impulse [s]',round(isp[0],2))

figure, axis = plt.subplots(3, 2)
  
axis[0, 0].plot(t, F)
axis[0, 0].set_title('Thrust [N] vs Time [s]')

  
axis[0, 1].plot(t, P/10**5)
axis[0, 1].set_title('Tank Pressure [bar] vs Time [s]')

  
axis[1, 0].plot(t, m_dot)
axis[1, 0].set_title('Mass Flow Rate [kg/s] vs Time [s]')

  
axis[1, 1].plot(t, T)
axis[1, 1].set_title('Temperature [K] vs Time [s]')

axis[2, 1].plot(t, isp)
axis[2, 1].set_title('Specific Impulse [s] vs Time [s]')

plt.show()



