import ajustador as aju
import numpy as np
from ajustador import drawing,nrd_fitness
import os

dirname='./'
#name of model xml file for optimization
model_set='Model_simple'
#name of experimental data, a simulation file in this case
exp_set='Model_simple'
#molecule to compare between 'experiments' and simulations
mol={'glu':['glu'], 'buf':['buf']}
#directory to store output during optimization
tmpdir='/tmp/fit'

# number of iterations, use 1 for testing
iterations=100
# default popsize=8, use 3 for testing
popsize=8

os.chdir(dirname)
exp = aju.xml.NeurordResult(exp_set)

#specify parameters to vary, either from ReactionScheme or InitialConditions
P = aju.xml.XMLParam
params = aju.optimize.ParamSet(P('buf_conc', 0, min=0, max=1000, xpath='//NanoMolarity[@specieID="buf"]'),
                               P('glu_rate', 0, min=0, max=1, xpath='//Reaction[@id="glu--glubuf_id"]/forwardRate'))

#this command indicates that experiments are from a previous simulation
###################### END CUSTOMIZATION #######################################

fitness = nrd_fitness.specie_concentration_fitness(species_list=mol)

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, sigma=0.3)

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
for i,p in enumerate(fit.params.unscale(fit.optimizer.result()[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(fit.optimizer.result()[6])[i])

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

