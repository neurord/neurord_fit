import ajustador as aju
import numpy as np
from ajustador import drawing,neurord_fit

#name of model xml file for optimization
model_set='Model_simple'
#name of experimental data, a simulation file in this case
exp_set='Model_simple'
#molecule to compare between 'experiments' and simulations
mol=['glu']
#directory to store output during optimization
tmpdir='/tmp/fit'

# number of iterations, use 1 for testing
iterations=2#100
# default popsize=8, use 3 for testing
popsize=3


#specify parameters to vary
P = aju.xml.XMLParam
params = aju.optimize.ParamSet(P('glu_fwd_rate', 0, min=0, max=1, xpath='//Reaction[@id="glu--glubuf_id"]/forwardRate'),
                               P('glu_rev_rate', 0, min=0, max=1, xpath='//Reaction[@id="glu--glubuf_id"]/reverseRate'))

#this command indicates that experiments are from a previous simulation
exp = aju.xml.NeurordResult(exp_set)
###################### END CUSTOMIZATION #######################################

fitness = neurord_fit.specie_concentration_fitness(species_list=mol)

fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(iterations, sigma=0.3)

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)
