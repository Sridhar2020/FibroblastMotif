#!/usr/bin/python

import os
import numpy as np
from os import chdir as cd

commandline = []

for rate in np.array([300,400,500,600]):
            for delay in np.array([0,10,25]): 
                         for strong in np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]):
                                       for weak in np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]):
                                                   folder = "Motif3_ShallowSlope_HomogenFibMyo_EqualT1T2_{0}ms_Delay{1}ms_LocCoup_{2}nS_LongCoup_{3}nS".format(rate,delay,strong,weak)  
                                                   commandline.append("./set.sh BCL1 {0}".format(rate))
                                                   commandline.append("./set.sh DELAY {0}".format(delay))
                                                   commandline.append("./set.sh G_STRONG1 {0}".format(strong))
                                                   commandline.append("./set.sh G_WEAK1 {0}".format(weak))
                                                   commandline.append("mkdir {0}".format(folder))
                                                   commandline.append("./makefile")
                                                  # commandline.append("./aout >> Output")
                                                  #commandline.append("mv *.stf {0}".format(folder))
                                                  #commandline.append("mv *.txt {0}".format(folder))
                                                  #commandline.append("mv Output {0}".format(folder))     
                                                   commandline.append("mv aout {0}".format(folder))
                                                  #commandline.append("cp Run_myjob.sh {0}".format(folder))
                                                  #commandline.append("cp checkpoint575000.out {0}".format(folder))
                                                  #commandline.append("cd {0}/".format(folder)) 
                                                  #commandline.append("qsub Run_myjob.sh")
                                                  #commandline.append("cd ..")            
                                                  #commandline.append("mv STFFiles/*.stf.gz {0}".format(folder)) 
                                                   
for command in commandline:
    print(command)
    os.system(command)
