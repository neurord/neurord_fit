#How to use :doc:`ajustador` to fit a NeuroRD model of CamKII activation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ajustador as aju
import numpy as np
from ajustador import drawing, loadconc, nrd_fitness
from ajustador.helpers import converge ,save_params
import os

dirname='camkii/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set='Model-CKnew4p-ss'
exp_set='Model-CKnew-ss0_1.5x' #set of data files corresponding to model files; files may contain several molecules
mol={"CaMKII": ['CKpCamCa4','CKCamCa4','CKp',]} #which molecule(s) to match in optimization
tmpdir='/tmp/PP1_4p'+dirname 
os.chdir(dirname)

#this command indicates that experiments are from a previous simulation
exp = aju.xml.NeurordResult(exp_set)

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=100
popsize=8
test_size=25 #for convergence

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('CKpCamPP1_fwd_rate', 4e-9, min=0, max=1e-3, xpath='//Reaction[@id="CKpCamPP1_bind"]/forwardRate'),
                               P('CKpCamPP1_bak_rate', 0.34e-3, min=0, max=1e-3, xpath='//Reaction[@id="CKpCamPP1_bind"]/reverseRate'),
                               P('CKpCamPP1_kcat_rate', 0.086e-3, min=0, max=1e-3, xpath='//Reaction[@id="CKpCamPP1_reac"]/forwardRate'))

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

#print centroid [0] and stdev [6] of cloud of good results:
for i,p in enumerate(fit.params.unscale(fit.optimizer.result()[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(fit.optimizer.result()[6])[i])

save_params.save_params(fit,0,1)

########################################## Next model
model_set='Model-CKnew-ss' #uses 6 parameter optimization
tmpdir='/tmp/PP1_6p'+dirname 

fit2 = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit2.load()
fit2.do_fit(iterations, popsize=popsize,sigma=0.3)
mean_dict2,std_dict2,CV2=converge.iterate_fit(fit2,test_size,popsize)

########################################### Done with fitting

#to look at fit history
aju.drawing.plot_history(fit2,fit2.measurement)
save_params.save_params(fit2,0,1)

#to look at centroid [0] or stdev [6] of cloud of good results:
for i,p in enumerate(fit2.params.unscale(fit2.optimizer.result()[0])):
    print(fit2.param_names()[i],'=',p, '+/-', fit2.params.unscale(fit2.optimizer.result()[6])[i])

'''
PP1 tuning
1. 4 parameters:
CKpCamPP1_fwd_rate = 1.39237415955e-11 +/- 1.25672092654e-11
CKpCamPP1_bak_rate = 0.000295587694727 +/- 3.50658584364e-06
CKpCamPP1_kcat_rate = 9.20116953403e-05 +/- 2.62340955889e-07

4p repeated - tune to Model-CKnew-ss0_1.5x.h5
CKpCamPP1_fwd_rate = 3.33693931971e-10 +/- 6.53396185018e-12
CKpCamPP1_bak_rate = 0.000597952563552 +/- 8.0341432889e-07
CKpCamPP1_kcat_rate = 8.96742999091e-05 +/- 4.14108096386e-08

6p tuned to previous best (flatter basal): Model-CKnew-ss0_1.5x.h5

CKpCamPP1_fwd_rate = 8.68282379153e-10 +/- 1.54763645414e-13
CKpCamPP1_bak_rate = 0.000560009040118 +/- 3.27005947467e-08
CKpCamPP1_kcat_rate = 8.7402349412e-05 +/- 1.08023999283e-09

Notice that Kb and Kcat quite similar between 4p and 6p
'''
