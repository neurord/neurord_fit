#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This demonstration fits a single reaction model (Model_simple.xml) to itself.
#to use this to optimize rate constants for a reaction, do the following steps
#1. create Model with the reaction (or reactions) to optimize
#2. create an array of value versus time with format the same as exp.output.counts().loc[0, :, 'glu', 0]
#3. update pop2= in the fitness definition
#4. update the pop1 - pop2 comparison to match time samples

import ajustador as aju
import tempfile
from ajustador import drawing

#edit the path as needed
model = aju.xml.open_model('neurord_fit/Model_simple.xml')
exp = aju.xml.NeurordResult('neurord_fit/Model_simple.h5')

P = aju.xml.XMLParam
params = aju.optimize.ParamSet(P('glu_fwd_rate', 0, min=0, max=1, xpath='//Reaction[@id="glu--glubuf_id"]/forwardRate'),P('glu_rev_rate', 0, min=0, max=1, xpath='//Reaction[@id="glu--glubuf_id"]/reverseRate'))

# # Just a test
# model2 = aju.xml.update_model(model, params)
# aju.xml.write_model(model2, '/dev/stdout')
# print()

# sim = aju.xml.NeurordSimulation('/tmp', model=model, params=params)

def specie_concentration_fitness(*, voxel=0, species, trial=0):
    def fitness(sim, measurement, full=False):
        pop1 = sim.output.counts().loc[voxel, :, species, trial]
        pop2 = measurement.output.counts().loc[voxel, :, species, trial]
        # Note: without interpolation this only works when the timepoints are identical
        diff = pop2 - pop1
        if full:
            return diff
        else:
            return float((diff**2).mean()**0.5)
    return fitness

fitness = specie_concentration_fitness(species='glu')
#print(fitness(sim, exp))

# first fit, just do a single iteraction, default popsize (=8)
fit = aju.optimize.Fit('/tmp/out1', exp, model, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
fit.do_fit(100, sigma=0.3)

#to look at parameters:
popsize=8
for i in range(popsize):
  print(fit[i].params)

  #to look at centroid [0] or stdev [6] of cloud of good results:
fit.params.unscale(fit.optimizer.result()[0])
fit.params.unscale(fit.optimizer.result()[6])

#to recall names of parameters
fit.param_names()

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#more examples
# smaller population (popsize=3)
fit2 = aju.optimize.Fit('/tmp/out2', exp, model, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit2.load()
#fit2.do_fit(200, sigma=0.3, popsize=3)

# popsize=8, no seed set in Model_simple.xml
fit3 = aju.optimize.Fit('/tmp/out3', exp, model, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit3.load()
fit3.do_fit(200, sigma=0.3, popsize=8)
