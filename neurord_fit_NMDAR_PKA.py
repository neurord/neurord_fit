#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing, loadconc, nrd_fitness
from ajustador.helpers import converge ,save_params
import os

dirname = 'nmdar_pka/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set = 'Model-NMDAR-PKA'
exp_set = 'pNMDAR_percent' #set of data files corresponding to model files; files may contain several molecules
mol = {'NMDAR': ['pNMDAR']} #which molecule(s) to match in optimization
tmpdir = '/tmp/NMDAR'
norm_method='percent'
os.chdir(dirname)

#this command indicates that experiments are from a previous simulation
exp = loadconc.CSV_conc_set(exp_set, stim_time=10)
# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations = 10
popsize = 3
test_size = 0 #for convergence

P = aju.xml.XMLParam
print(P)
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('pNMDARPP1_fwd_rate',
                                 4e-9, min=0, max=1e-3,
                                 xpath='//Reaction[@id="bindpNMDARPP1"]/forwardRate'),
                               P('pNMDARPP1_back_rate',
                                 0.34e-3, min=0, max=1e-3,
                                 xpath='//Reaction[@id="bindpNMDARPP1"]/reverseRate'),
                               P('pNMDARPP1_kcat_rate',
                                 0.086e-3, min=0, max=1e-3,
                                 xpath='//Reaction[@id="reacpNMDARPP1"]/forwardRate'))

###################### END CUSTOMIZATION #######################################

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol,
                                                   norm=norm_method)
fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
print(fit.model)
fit.do_fit(iterations, popsize=popsize, sigma=0.3)
#fit.do_fit(iterations, popsize=popsize, seed=62839)
#mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting, look at results

#to look at fit history
aju.drawing.plot_history(fit, fit.measurement, Norm='percent')

#to look at centroid [0] or stdev [6] of cloud of good results:
if callable(fit.optimizer.result):
        result = fit.optimizer.result()
    else:
        result = fit.optimizer.result
for i,p in enumerate(fit.params.unscale(result[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(result[6])[i])

save_params.save_params(fit,0,1)
