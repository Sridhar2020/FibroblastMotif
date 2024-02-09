#!/usr/bin/python

import os
import numpy as np

commandline = []

for rate in np.array([300,400,500,600]):
            for delay in np.array([0,10,25]): 
                         for strong in np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]):
                                       for weak in np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]):
                                                   folder = "Motif3_SteepSlope_HomogenFibMyo_EqualT1T2_{0}ms_Delay{1}ms_LocCoup_{2}nS_LongCoup_{3}nS".format(rate,delay,strong,weak)  
                                                   os.system('cp Run_myjob.sh {0}'.format(folder)) 
                                                   os.chdir(folder)
                                                   #print(os.getcwd())
                                                   #os.system('qsub Run_myjob.sh')
                                                   os.system('./aout >> Output')
                                                   os.chdir('../')
                                                   #print(os.getcwd())
                                                   
for command in commandline:
    print(command)
    os.system(command)
