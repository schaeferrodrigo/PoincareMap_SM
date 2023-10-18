import parameters as par
from basic_functions.functions import *
from basic_functions.tau_star import *
#from integ_mthd.ReturnMap import *
from basic_functions.vector_field import *
from scipy import integrate
import os

def main():

    dir = os.getcwd()

    command = 'export PYTHONPATH="' +dir+'"'
    print(command)
    os.system(command)

    #values of mu parameters
    a = par.a

    I1_0 , I2_0  = par.I_1_initial ,  par.I_2_initial

    domainTheta2 = par.domainTheta2 # interval of values for theta_2

    #Energy level
    # note that for I =(0,0). the value tau is always tau = 0
    energyLevel = RPF(a, par.I_energy , par.thetaEnergy1, par.thetaEnergy2, 0)
    print('fixed energy level = ' , energyLevel)

    maxSteps = par.maxSteps
    maxInter = par.maxInter

    file = open("data/energyConservation", "w")


    for theta2_0 in domainTheta2:

        nameFile = str(theta2_0)

        errorEnergy = 0
        errorCrest = 0
        differenceEnergy= 0.
        differenceCrest= 0.

        PMfile = open("dataPM/"+ nameFile, "w")


        theta1_0 =  theta_1_init(energyLevel, a, theta2_0)

        print("initial angles are theta_1 = ", theta1_0, 'and theta_2 = ', theta2_0 )

        iniCond = np.array([I1_0, I2_0,theta1_0, theta2_0 , 0])


        solution = integrate.DOP853(F , 0. , iniCond, t_bound = maxSteps, atol = 1e-12)

        integrationTime = np.array([0.])

        step = 1
        numberInter = 0
        oldI1,oldI2 , oldTheta1,oldTheta2 = I1_0, I2_0,theta1_0, theta2_0

        while  step < maxSteps and numberInter < maxInter:
        #while  step <= maxSteps and I_1_list[-1]*I_1_list[-2] >=0:

            solution.step()

            integrationTime = np.append(integrationTime , solution.t)

            I_1 , I_2, theta_1, theta_2, tau = solution.y

            print('I_1 = ', I_1)
            print('step =', step)
            print('numberInter = ', numberInter)

            if I_1*oldI1 <0 and oldI1 <0:
                if np.abs(I_1) < 1e-12 or np.abs(oldI1) <1e-12:
                    if np.abs(I_1) < 1e-12:
                        Sigma_I_2 , Sigma_theta_2 = I_2 , theta_2
                    else:
                        Sigma_I_2 , Sigma_theta_2 = oldI2 , oldTheta2
                else:
                    t_i , t_f = integrationTime[-2] , integrationTime[-1]
                    I_i , I_f = oldI1 , I_1
                    I_c = I_1
                    sol = np.zeros(4)
                    while np.abs(I_c) > 1e-12:
                        t_c = t_i + (t_f-t_i)/2.
                        sol = solution.dense_output().__call__(t_c)
                        I_c = sol[0]
                        if I_i*I_c <0:
                            I_f , t_f = I_c , t_c
                            #print('t_f =', t_f)
                        else:
                            I_i , t_i = I_c , t_c
                            #print('t_i= ', t_i)
                        Sigma_I_2, Sigma_theta_2 = sol[1] , sol[3]
                    #print('I_c = ',I_c)

                PMfile.write(str(Sigma_theta_2)+' '+str(Sigma_I_2)+'\n')
                step =1
                numberInter += 1
                iniCond = np.array([I_1 , I_2, theta_1, theta_2, tau])
                solution = integrate.DOP853(F , 0. , iniCond, t_bound = maxSteps,atol = 1e-12)
            else:
                step +=1
            oldI1,oldI2 , oldTheta1,oldTheta2 =I_1 , I_2, theta_1, theta_2

            if np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))> 1e-8 and errorEnergy == 0:
                errorEnergy = 1
                differenceEnergy = np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
                #print(differenceEnergy)

            if np.abs(crests(a, [I_1, I_2], theta_1, theta_2, tau)) > 1e-8 and errorCrest == 0:
                errorCrest = 1
                differenceCrest =  np.abs(crests(a, [I_1, I_2], theta_1, theta_2, tau))
                #print(differenceCrest)

            integrationTime = integrationTime[-2:]
        PMfile.close()


            # print('energy is conserved? ', energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
            # print('equations of crest is hold?', crests(a, [I_1, I_2], theta_1, theta_2,  tau))
            # print(solution.y)
        if errorEnergy == 1:
            file.write('energy is not conserved = ' + str(differenceEnergy)+'\n')
        if errorCrest == 1:
            file.write('crest is not satisfied = '+str(differenceCrest)+'\n')

    file.close()



if __name__ == '__main__':
    main()
