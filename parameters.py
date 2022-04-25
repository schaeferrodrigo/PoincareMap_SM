import numpy as np

#parameters to be used


a = [0.3 ,0.3] # the sum of the entries has to be less than 0.625

I_1_initial = 0.00001 # transversal section

I_2_initial = 0.0001 # fixed initial value

#domain of initial values of theta_2
delta = 0.3
domainTheta2 = np.linspace(np.pi + delta , 3.12, 1000)

#Fixed energy level
I_energy = [0,0]
theta_energy = 5*np.pi/4
#

#value of epsilon (it can be view as a step in time)
eps = 0.001
