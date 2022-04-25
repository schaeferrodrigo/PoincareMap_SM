import parameters as par
from basic_functions.functions import *
from basic_functions.tau_star import *

def main():
    #valeus of mu parameters
    a = par.a

    trans_section = par.I_1_initial # transversal section
    init_value_I_2 = par.I_2_initial
    domainTheta2 = par.domainTheta2 # initival value for theta_2

    #Energy level
    # note that for I =(0,0). the value tau is always tau = 0
    enr_lvl = RPF(a, par.I_energy , par.theta_energy, par.theta_energy, 0)
    print('fixed energy level = ' , enr_lvl)

    #time step
    eps = par.eps

    for theta_2_0 in domainTheta2:
        print((enr_lvl-A(1,1) - A(0,a[1])*np.cos(theta_2_0))/A(0,a[0]))
        theta_1_0 = np.arccos((enr_lvl-A(1,1) - A(0,a[1])*np.cos(theta_2_0))/A(0,a[0]))

        print(theta_1_0)



if __name__ == '__main__':
    main()
