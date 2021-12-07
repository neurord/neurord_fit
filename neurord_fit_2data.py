#How to use :doc:`ajustador` to fit a NeuroRD model to experimental data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a single reaction model (Model_pSynGap) to data.
'''

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge
import os

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='syngap_ras/'  #where data and model are stored.  Multiple datafiles allowed
model_set='Model_syngap_ras'
exp_name='walkup_JBC_2015' #name of data file selected from dirname; each file may contain several molecules
mol=['pSynGap','RasGDP'] #which molecule(s) to match in optimization
tmpdir='/tmp/'+dirname

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=25
popsize=8
test_size=25

os.chdir(dirname)
exp=loadconc.CSV_conc_set(exp_name)

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

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol)

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
mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
if callable(fit.optimizer.result):
    result=fit.optimizer.result()
else:
    result=fit.optimizer.result
for i,p in enumerate(fit.params.unscale(result[0])):
        print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(result[6])[i])

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

