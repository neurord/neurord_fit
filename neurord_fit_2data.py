#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a single reaction model (Model_pSynGap) to data.
Next steps:
1. fix plotting the traces in drawing.py - test
2. Create data set with two conditions - test 
3. Move fitness function to ajustador.fitnesses - test
'''

import ajustador as aju
import numpy as np
import loadconc
from ajustador import drawing

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='./walkup/'  #where is data stored.  Multiple datafiles allowed
model_set=dirname+'Model_syngap_ras'
exp_name=dirname+'walkup_JBC_2015' #name of data file selected from dirname; each file may contain several molecules
mol=['pSynGap','RasGDP'] #which molecule(s) to match in optimization
tmpdir='/tmp/'+dirname

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=1
popsize=3

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

###################### fitness_functions - move to nrd_fitness.py in ajustador #######################################
def specie_concentration_fitness(*, voxel=0, species_list, trial=0):
    def fitness(sim, measurement, full=False):
        fitarray=np.zeros((len(species_list),len(sim.output)))
        fit_dict={}
        for i,species in enumerate(species_list):
            fit_dict[species]={}
            for j,stim_set in enumerate(sim.output):
                pop1=aju.nrd_output.nrd_output_conc(stim_set,species)
                stim_set.__exit__()
                #pop1 = stim_set.concentrations().loc[voxel, :, species, trial]
                if isinstance(measurement,aju.xml.NeurordResult):
                    #pop2 = measurement.output[i].concentrations().loc[voxel, :, species, trial]
                    pop2 = aju.nrd_output.nrd_output_conc(measurement.output[j],species)
                    diff = pop2 - pop1
                    max_mol=np.mean([np.max(pop1.values),np.max(pop2.values)])
                else:  #measurement is experimental data, stored as CSV_conc_set
                    #print(measurement.data[j].name, type(measurement.data[j]))
                    pop2 = measurement.data[j].waves[species].wave
                    wave1y=pop1.values[:,0]
                    wave1x=pop1.index
                    # Note: np.interp(x1,x2,y2) returns values for y2 corresponding to x1 timepoints
                    pop1y=np.interp(pop2.x,wave1x,wave1y)
                    diff = pop2.y - pop1y
                    max_mol=np.mean([np.max(pop1.values),np.max(pop2.y)])
                diffnorm = diff if max_mol==0 else diff/max_mol
                fit_dict[species][stim_set.injection]=float((diffnorm**2).mean()**0.5)
                fitarray[i][j]=float((diffnorm**2).mean()**0.5)
        fitness=np.mean(fitarray)
        #print ('fitarray', fitarray)
        if full:
            return fit_dict
        else:
            return fitness
    return fitness

fitness = specie_concentration_fitness(species_list=mol)

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

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))


#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

'''
#to plot data (clicking on point in fit history doesn't work because no "injection")
#if other than last sim is desired, assign that value to fitnum
def plot_traces(fitX,mollist,fitnum=-1):
    from matplotlib import pyplot as plt
    import os
    plt.ion()
    bestdir=fitX[fitnum].tmpdir.name
    fname=bestdir+'/model.h5'
    print ('PLOTTING', bestdir, 'parameters: ',fit[fitnum])
    fig,axes=plt.subplots(len(mollist),1,sharex=True,figsize=(6,9))
    fig.canvas.set_window_title(fitX.name)
    best=aju.nrd_output.Output(fname)
    for i,mol in enumerate(mollist):
        simdata=nrd_output_conc(best,mol)
        axes[i].plot(simdata.index,simdata.values[:,0])
        axes[i].plot(exp.waves[mol].wave.x,exp.waves[mol].wave.y)
        axes[i].set_ylabel(mol)
    axes[i].set_xlabel('time, msec')

fitnum=-1
plot_traces(fit,mol)
'''
