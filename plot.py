import numpy as np
import matplotlib.pyplot as plt
import os
import parameters as par

def main():

    files = os.listdir('dataPM/')
    plt.figure()
    plt.axis([0,2*np.pi , -0.04,0.04])
    for file in files:

        theta2_values = np.array([])
        I2_values = np.array([])

        arx =  open('dataPM/'+file,'r')

        lines = arx.readlines()
        if len(lines) != 0:
            for line in lines:
                theta2,I2 = line.split()
                theta2,I2 = float(theta2)%2*np.pi,float(I2)
                theta2_values = np.append(theta2_values, theta2)
                I2_values = np.append(I2_values,I2)
            plt.plot(theta2_values,I2_values,'.', markersize= .5)

        arx.close()
    mus =str(int(par.a[0]*100))+'_'+ str(int(par.a[1]*100))
    plt.xlabel(r'$\theta_2$' , size = 15)
    plt.ylabel(r'$I_2$' , rotation = 0, size = 15 )
    print(mus)
    name_fig = 'figs/pm_'+mus
    print(name_fig)
    plt.savefig(name_fig)
    #plt.show()


if __name__ == '__main__':
    main()
