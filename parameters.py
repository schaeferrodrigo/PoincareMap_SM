import numpy as np

#parameters to be used


a = [0.2,0.3] # the sum of the entries has to be less than 0.625

I_1_initial = 1e-17 # transversal section

I_2_initial = 1e-17  # fixed initial value

#Fixed energy level
I_energy = [0,0]
thetaEnergy1 = 5*np.pi/4
thetaEnergy2 = 5*np.pi/4
#

#domain of initial values of theta_2
delta = 0.2
#domainTheta2 = np.linspace(3.9,4, 10)
numberPoints = 500
domainTheta2 = np.linspace(thetaEnergy2-delta, thetaEnergy2 + delta, numberPoints)

#integration method
maxSteps = 1e6

# number of intersections of tranversal sections:
maxInter = 1000
