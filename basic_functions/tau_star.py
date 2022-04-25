from basic_functions.functions import *

import numpy as np

def bissec_method( a ,I, theta, tau_1 , tau_2  , tol = 1e-2 ):
    theta_1, theta_2 = theta
    while np.abs( tau_1 - tau_2 ) > tol:
        tau_c = ( tau_1 + tau_2 )/2
        if crests( a ,I, theta_1, theta_2,  tau_1 ) * crests(a,I,theta_1 , theta_2,  tau_c) < 0:
            tau_2 = tau_c
        else:
            tau_1 = tau_c
    return [tau_1 , tau_2]

def secant_method( a, I,theta, tau_1 , tau_2 ,tol = 1e-10):
    theta_1,theta_2 = theta
    tau_sec = tau_1 - crests(a ,I , theta_1 , theta_2,  tau_1 )*( tau_1 - tau_2)/(crests( a ,I , theta_1 , theta_2,  tau_1) - crests( a , I , theta_1 , theta_2,  tau_2) )
    while np.abs( crests(a , I, theta_1 , theta_2,  tau_sec) ) > tol:
        if crests( a ,I , theta_1 , theta_2,  tau_1) * crests(a, I , theta_1 , theta_2,  tau_sec) < 0 :
            tau_2  = tau_sec
        else:
            tau_1 = tau_sec
        tau_sec =  tau_1 - crests(a ,I , theta_1 , theta_2,  tau_1)*( tau_1 - tau_2)/(crests(a , I , theta_1 , theta_2,  tau_1 ) -crests( a ,I , theta_1 , theta_2,  tau_2 ) )
    return tau_sec

def tau(a ,I ,theta ,sign , tau_initial = 0 , step = 0.05):
    theta_1,theta_2 = theta
    tau_2 = tau_initial
    if sign == 'neg':
        tau_1 = -step
    else:
        tau_1 = step
    sign = np.sign(tau_1)
    while crests(a , I , theta_1 , theta_2,  tau_1 ) * crests(a ,I , theta_1 , theta_2,  tau_2 ) > 0 :
        tau_2 = tau_1
        tau_1 = tau_1 + sign * step
    candidate = bissec_method(a , I,theta, tau_1 , tau_2  )
    tau_1 = candidate[0]
    tau_2 = candidate[1]
    tau_star = secant_method( a , I,theta, tau_1 , tau_2) #método secanda
    return tau_star

def assign_tau(a, I, theta ):
    value_of_tau_pos = tau(a, I, theta , 'pos' )
    value_of_tau_neg = tau( a,I , theta , 'neg' )
    if np.minimum( np.abs(value_of_tau_neg) , np.abs(value_of_tau_pos) ) == np.abs(value_of_tau_neg):
        return value_of_tau_neg
    else:
        return value_of_tau_pos
