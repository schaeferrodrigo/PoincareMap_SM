import  numpy as np


def alpha(I):
    if I != 0:
        return (I**2)*np.sinh(np.pi/2)/np.sinh(I*np.pi/2)
    else:
        return 0


def A(I ,a ):
    if I != 0 :
        return 2*np.pi* I *a/(np.sinh(np.pi *I/2))
    else:
        return 4*a

def der_A( I , theta):
    return np.pi*a*(2 - np.pi* I*(1/np.tanh(np.pi*I/2))) * (1/np.sinh(np.pi*I/2))

def RPF(a ,I , theta_1 ,theta_2, tau):
    l= A(I[0] , a[0]) * np.cos(theta_1-I[0]*tau) + A(I[1], a[1])* np.cos(theta_2 - I[1]*tau) +A(1,1)*np.cos(-tau)
    return l

def crests(a, I, theta_1, theta_2,  tau):
    eq = alpha(I[0])*a[0]*np.sin(theta_1 - I[0]*tau)+ alpha(I[1])*a[1]*np.sin(theta_2-I[1]*tau) + np.sin(-tau)
    return eq

def der_I_RPF(a,I, theta , tau):
    return -A(I)*np.sin(theta - I*tau)

def der_theta_RPF(a,I,theta, tau):
    return -(der_A(I, theta)*np.cos(theta-I*tau) + tau*A(I,a)*np.sin(theta - I*tau))
