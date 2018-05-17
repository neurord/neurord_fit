#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a model of CamKII to a previous models results
   simulate CamCa4 pulses at different frequencies, possibly include models with 6 s of constant stimulation at different CamCa4 concentrations.  May need to inject Ca and CaBuf instead of CamBuf to get better CamCa4 levels. 
Next Steps:
4. Add variation of IC parameters (important)

5. Specify parameters to vary in a separate file (less important)
     possibly with xml specifications similar to the model instead of xpath - see neurord3_params_sims.py and neurord_2params_sims.py

6. add features to fitness function, e.g. peak amplitude and width 
    need to specify stimulus onset time, as with injection current. 
    Attach fitness values to fit object?
'''

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,neurord_fit
from ajustador.helpers import converge

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='camkii/'  #location of data and model files.
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set=dirname+'Model-CKnew-Cahz'
exp_set=dirname+'Model-CKold-Cahz' #set of data files corresponding to model files; files may contain several molecules
mol=['CKpCamCa4','CKCamCa4'] #which molecule(s) to match in optimization
tmpdir='/tmp/'+dirname 

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=100
popsize=8
test_size=25 #for convergence

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('CK2_fwd_rate', 0, min=0, max=1e-6, xpath='//Reaction[@id="CKCam_pow2"]/forwardRate'),
                               P('CK4_fwd_rate', 0, min=0, max=1e-9, xpath='//Reaction[@id="CKCam_pow4"]/forwardRate'),
                               P('CK3_fwd_rate', 0, min=0, max=1e-12, xpath='//Reaction[@id="CKCam_pow3"]/forwardRate'),
                               P('CK1_CKp2_fwd_rate', 0, min=0, max=1e-6, xpath='//Reaction[@id="CK1_CKpCam_pow2"]/forwardRate'),
                               P('CK2_CKp1_fwd_rate', 0, min=0, max=1e-6, xpath='//Reaction[@id="CK2_CKpCam_pow1"]/forwardRate'),
                               P('CK2_CKp2_fwd_rate', 0, min=0, max=1e-9, xpath='//Reaction[@id="CK2_CKpCam_pow2"]/forwardRate'))

###################### END CUSTOMIZATION #######################################
exp = aju.xml.NeurordResult(exp_set)

fitness = neurord_fit.specie_concentration_fitness(species_list=mol)

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

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)
