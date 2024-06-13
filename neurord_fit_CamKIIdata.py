#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os

dirname='camkii/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='Model-CKnew-Cahz'
exp_set='CamKII_Hz' #set of data files corresponding to model files; files may contain several molecules
mol={"CKpCamCa4": ['CKpCamCa4']} #which molecule(s) to match in optimization
tmpdir='/tmp/st0_4params'+dirname 
os.chdir(dirname)

# Use loadconc.CSV_conc_set if data to match are csv format (typically from wet experiments)
exp = loadconc.CSV_conc_set(exp_set)

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=100
popsize=8
test_size=25 #for convergence

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('CK2_fwd_rate', 1e-9, min=0, max=1e-6, xpath='//Reaction[@id="CKCam_pow2"]/forwardRate'),
                               P('CK4_fwd_rate', 1e-15, min=0, max=1e-12, xpath='//Reaction[@id="CKCam_pow4"]/forwardRate'),
                               P('CK3_fwd_rate', 1e-12, min=0, max=1e-9, xpath='//Reaction[@id="CKCam_pow3"]/forwardRate'),
                               P('CK2_CKp2_fwd_rate', 1e-12, min=0, max=1e-9, xpath='//Reaction[@id="CK2_CKpCam_pow2"]/forwardRate'))

###################### END CUSTOMIZATION #######################################

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol)

############ Test fitness function
#model=dirname+'Model-CKnew-Cahz1.xml'
#sim = aju.xml.NeurordSimulation('/tmp', model=model, params=params)
#sim2=aju.xml.NeurordResult('Model_syngap_ras.h5')
#print(fitness(sim2, exp))
################

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, popsize=popsize,sigma=0.3)
mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)
if callable(fit.optimizer.result):
    result = fit.optimizer.result()
else:
    result = fit.optimizer.result
#print centroid [0] and stdev [6] of cloud of good results:
for i,p in enumerate(fit.params.unscale(result[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(result[6])[i])

save_params.save_params(fit,0,1)

'''
a. try starting from previous good parameters  - this worked!          previous:                           new:old
CK2_fwd_rate = 3.43460459241e-10 +/- 7.36199078809e-12                 CK2_fwd_rate', 4.78e-10 - similar - 71%
CK4_fwd_rate = 2.96293987538e-16 +/- 5.22679398079e-18                 CK4_fwd_rate', 2.40e-16 - similar - 123%
CK3_fwd_rate = 6.38652104333e-15 +/- 1.56822600609e-15                 CK3_fwd_rate', 1.38e-14 - diff  -   46%
CK1_CKp2_fwd_rate = 2.67337675692e-13 +/- 1.92768856981e-14            CK1_CKp2_fwd_rate', 1.79e-13 - similar - 149%
CK2_CKp1_fwd_rate = 1.32611014657e-14 +/- 4.58697020871e-14            CK2_CKp1_fwd_rate', 1.10e-12 - HUGH diff - 1%, large var
CK2_CKp2_fwd_rate = 1.75527421788e-17 +/- 6.5301636154e-19             CK2_CKp2_fwd_rate', 2.20e-17 - similar - 79%

b. opt to all data, starting from 0. DOESN"T WORK
c. opt to data, without 6s data, starting from 0: DOESN"T WORK
d. opt to data, with 6s data, but start from coarse values - minimum fitness of 0.1
CK2_fwd_rate = 3.83119676007e-10 +/- 9.8597395351e-13
CK4_fwd_rate = 2.23694832822e-16 +/- 1.54061581795e-18
CK3_fwd_rate = 3.55792141502e-13 +/- 2.20969578794e-15
CK1_CKp2_fwd_rate = 3.0330636368e-13 +/- 2.29809824284e-13
CK2_CKp1_fwd_rate = 2.39399544094e-13 +/- 2.71561338665e-13 - STILL LARGE VARIANCE
CK2_CKp2_fwd_rate = 1.10450093247e-18 +/- 1.41461185682e-16 - EVEN LARGER VARIANCE (but small in previous sim)

Next: opt to data, exclude CK2_CKp1_fwd_rate (set = 0) - similar minimum fitness of 0.1
CK2_fwd_rate = 1.25915358328e-09 +/- 2.55635426689e-12
CK4_fwd_rate = 1.55853698323e-16 +/- 4.27641783607e-18
CK3_fwd_rate = 3.7001965304e-13 +/- 5.81743571691e-15
CK1_CKp2_fwd_rate = 5.60003015714e-13 +/- 4.76116161324e-13 - variance as large as value
CK2_CKp2_fwd_rate = 1.16831852933e-17 +/- 4.56376338119e-16

Next: opt to data, exclude CK1_CKp2_fwd_rate (set = 0)  - 4 parameters
CK2_fwd_rate = 7.89325641207e-10 +/- 1.20234737985e-12
CK4_fwd_rate = 7.50031075336e-17 +/- 1.92187649653e-19
CK3_fwd_rate = 9.66232673723e-13 +/- 3.90664628821e-16
CK2_CKp2_fwd_rate = 1.07481418709e-16 +/- 1.15866464049e-17

Next, using 4 and 6 param results, repeat PP1 optimization (update Rxn_CamKIInew_Ca.xml)

'''
