#How to use :doc:`ajustador` to fit a NeuroRD model to experimental data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a single reaction model (Model_pSynGap) to data.
'''

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,neurord_fit
from ajustador.helpers import converge

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='syngap_ras/'  #where is data stored.  Multiple datafiles allowed
model_set=dirname+'Model_syngap_ras'
exp_name=dirname+'walkup_JBC_2015' #name of data file selected from dirname; each file may contain several molecules
mol=['pSynGap','RasGDP'] #which molecule(s) to match in optimization
tmpdir='/tmp/'+dirname

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=1
popsize=3
test_size=25

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('phos_fwd_rate', 0, min=0, max=1, xpath='//Reaction[@id="CKpCamCa4+SynGap--CKpCamCa4SynGap"]/forwardRate'),
                               P('phos_rev_rate', 0, min=0, max=1, xpath='//Reaction[@id="CKpCamCa4+SynGap--CKpCamCa4SynGap"]/reverseRate'),
                               P('phos_kcat_rate', 0, min=0, max=1, xpath='//Reaction[@id="CKpCamCa4SynGap--CKpCamCa4+pSynGap"]/forwardRate'),
                               P('gap_kf_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTP+SynGap--RasGTPGap"]/forwardRate'),
                               P('gap_kb_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTP+SynGap--RasGTPGap"]/reverseRate'),
                               P('gap_kcat_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTPGap--SynGap+RasGDP"]/forwardRate'),
                               P('Pgap_kf_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTP+pSynGap--RasGTPGap"]/forwardRate'),
                               P('Pgap_kb_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTP+pSynGap--RasGTPGap"]/reverseRate'),
                               P('Pgap_kcat_rate', 0, min=0,max=1,xpath='//Reaction[@id="RasGTPGap--pSynGap+RasGDP"]/forwardRate'))

###################### END CUSTOMIZATION #######################################

exp=loadconc.CSV_conc_set(exp_name)

fitness = neurord_fit.specie_concentration_fitness(species_list=mol)

############ Test fitness function
#sim = aju.xml.NeurordSimulation('/tmp', model=model, params=params)
#cp /tmp/???/model.h5 modelname.split('.')[0]+'.h5'
#sim2=aju.xml.NeurordResult('Model_syngap_ras.h5')
#print(fitness(sim2, exp))
################

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, popsize=popsize,sigma=0.3)
#mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))


#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

