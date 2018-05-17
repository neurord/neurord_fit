#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a NE to cAMP to Epac_FRET (with PKA&PDE) to data.
Note that stimulation must occur at same time in simulation and experiment
Necessary to specify time of stimulation, in seconds
Eliminate calcium & calmodulin reactions to speed simulation, and decrease num comps in 6.5u from 9 to 5.
'''

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,neurord_fit
from ajustador.helpers import converge
import os

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='fret_cAMP/'  #where data and model file are stored.  Multiple datafiles allowed
model_set='Model_VinceCaPKAsubset_ACsame-dia'
exp_name='FretPercent' #name of data file selected from dirname; each file may contain several molecules
mol=['Epac1cAMP'] #which molecule(s) to match in optimization
tmpdir='/tmp/fixed_'+dirname
start_stim=100  #time of onset of stimulation, in seconds; should be able to extract from model files; stim time must match data
norm_method='percent' #convert molecule concentration into a percent change from baseline.  Use for comparing to FRET data

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=50
popsize=6 # reduce from 8 to avoid saturating processors for longer neurord simulations
test_size=25

P = aju.xml.XMLParam
#list of parameters to change/optimize
#improvement: allow rate (multiply all kf,kb,kcat) or KM (vary what? Kb?)
params = aju.optimize.ParamSet(P('phospLRGs_fwd_rate', 0, min=0, max=1, xpath='//Reaction[@id="phospLRGs1"]/forwardRate'),
                               P('phospLRGs_bak_rate', 0, min=0,max=1,xpath='//Reaction[@id="phospLRGs1"]/reverseRate'),
                               P('phospLRGs_cat_rate', 0, min=0,max=1,xpath='//Reaction[@id="phospLRGs2"]/forwardRate'),
                               P('dphospD1R_fwd_rate', 0, min=0,max=1,xpath='//Reaction[@id="dphospD1R"]/forwardRate'),
                               P('GasGTP_hydrolysis', 0, min=0,max=1,xpath='//Reaction[@id="GasGTP_disso"]/forwardRate'),
                               P('AC1_GasGTP_GAP', 0, min=0,max=1,xpath='//Reaction[@id="AC1_GasGTP_GAP"]/forwardRate'),
                               P('dphospPDE4_fwd_rate', 0, min=0,max=1,xpath='//Reaction[@id="dphospPDE4"]/forwardRate'))

###################### END CUSTOMIZATION #######################################
os.chdir(dirname)
exp=loadconc.CSV_conc_set(exp_name)
fitness = neurord_fit.specie_concentration_fitness(species_list=mol,start=start_stim,norm=norm_method)
fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, popsize=popsize,sigma=0.3)
#mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting, look at results

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))


#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

