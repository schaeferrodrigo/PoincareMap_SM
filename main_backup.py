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


    I_1_list = np.array([0.])
    I_2_list = np.array([0.])
    theta_1_list = np.array([0.])
    theta_2_list =np.array([0.])

    #values of mu parameters
    a = par.a

    I_1_list =  np.append(I_1_list,par.I_1_initial) # transversal section
    I_2_list = np.append(I_2_list,par.I_2_initial)
    domainTheta2 = par.domainTheta2 # intereval of values for theta_2

    #Energy level
    # note that for I =(0,0). the value tau is always tau = 0
    energyLevel = RPF(a, par.I_energy , par.thetaEnergy1, par.thetaEnergy2, 0)
    print('fixed energy level = ' , energyLevel)

    PMfile = open("dataPM/pointsPoincare", "w")

    for theta_2_0 in domainTheta2:
        nameFile = str(theta_2_0)

        errorEnergy = 0
        errorCrest = 0
        differenceEnergy= 0.
        differenceCrest= 0.


        file = open("data/"+ nameFile, "w")



        theta_2_list = np.append(theta_2_list,theta_2_0)
        theta_1_list = np.append(theta_1_list, theta_1_init(energyLevel, a, theta_2_0))

        print("initial angles are theta_1 = ", theta_1_list[1], 'and theta_2 = ', theta_2_0 )

        iniCond = np.array([I_1_list[1], I_2_list[1],theta_1_list[1], theta_2_list[1] , 0])

        maxSteps = par.maxSteps
        solution = integrate.DOP853(F , 0. , iniCond, maxSteps ,rtol = 1e-2,atol = 1e-12)

        integrationTime = np.array([0.])

        step = 1
        while  I_1_list[-1]*I_1_list[-2] >=0:
        #while  step <= maxSteps and I_1_list[-1]*I_1_list[-2] >=0:
            solution.step()

            integrationTime = np.append(integrationTime , solution.t)

            I_1 , I_2, theta_1, theta_2, tau = solution.y
            I_1_list = np.append(I_1_list, I_1)
            I_2_list = np.append(I_2_list, I_2)
            theta_2_list = np.append(theta_2_list, theta_2)
            print('I_1 = ', I_1)
            #print('theta_2= ',   theta_2,' I_2 = ', I_2)
            print('step =', step)
            tauStar = assign_tau(a, [I_1,I_2], [theta_1,theta_2])
            #if np.abs(tau_star- tau)>1e-12:
            #    print('difference in tau =', np.abs(tau_star- tau))
            file.write(str(solution.y)+"\n")

            # crest = np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
            # starCrest = np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tauStar))
            # if np.abs(crest) <= np.abs(starCrest):
            #     print('tau')
            # else:
            #     print('tau_star')
            #     print(crest )
            #     print(starCrest)


            if np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))> 1e-8 and errorEnergy == 0:
                errorEnergy = 1
                differenceEnergy = np.abs(energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
                print(differenceEnergy)

            if np.abs(crests(a, [I_1, I_2], theta_1, theta_2, tau)) > 1e-8 and errorCrest == 0:
                errorCrest = 1
                differenceCrest =  np.abs(crests(a, [I_1, I_2], theta_1, theta_2, tau))
                print(differenceCrest)

            step +=1
            integrationTime = integrationTime[-2:]

            # print('energy is conserved? ', energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
            # print('equations of crest is hold?', crests(a, [I_1, I_2], theta_1, theta_2,  tau))
            # print(solution.y)
        if errorEnergy == 1:
            file.write('energy is not conserved = ' + str(differenceEnergy)+'\n')
        if errorCrest == 1:
            file.write('crest is not satisfied = '+str(differenceCrest)+'\n')

        file.close()


        I_i , I_f = I_1_list[-2] , I_1_list[-1]
        t_i , t_f = integrationTime[-2] , integrationTime[-1]

        if np.abs(I_i) < 1e-12 or np.abs(I_f) <1e-12:
            if np.abs(I_i) < 1e-12:
                I_2 , theta_2 = I_2_list[-2] , theta_2_list[-2]
            else:
                I_2 , theta_2 = I_2_list[-1] , theta_2_list[-1]
        else:
            I_c = I_1
            sol = np.zeros(4)
            while np.abs(I_c) > 1e-12:
                t_c = t_i + (t_f-t_i)/2.
                #print(t_f, t_i)
                #print('t_c = ',(t_f-t_i)/2.)
                #print(solution.dense_output().__call__(t_c))
                sol = solution.dense_output().__call__(t_c)
                I_c = sol[0]
                if I_i*I_c <0:
                    I_f , t_f = I_c , t_c
                    #print('t_f =', t_f)
                else:
                    I_i , t_i = I_c , t_c
                    #print('t_i= ', t_i)
                I_2, theta_2 = sol[1] , sol[3]
            print('I_c = ',I_c)

        PMfile.write(str(theta_2)+' '+str(I_2)+'\n')

        theta_2_list = theta_2_list[:1]
        theta_1_list = theta_1_list[:1]
        I_1_list = I_1_list[:2]
        I_2_list = I_2_list[:2]
        #print(I_1_list)

    PMfile.close()


        #returnPoint = ReturnM(a, iniCond, eps, 1000  )

        # I_1,I_2,theta_1, theta_2 = returnPoint
        # tau = assign_tau(a, [I_1, I_2], [theta_1, theta_2])
        #
        #
        # print('I_1 = ', I_1)
        # print('I_2 = ', I_2)
        # print('theta_1 = ', theta_1)
        # print('theta_2 = ', theta_2)
        # print('energy is conserved? ', energyLevel - RPF(a, [I_1 , I_2], theta_1, theta_2, tau))
        # print('equations of crest is hold?', crests(a, [I_1, I_2], theta_1, theta_2,  tau))


if __name__ == '__main__':
    main()
