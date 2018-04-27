#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a model of CamKII to a previous models results
   simulate CamCa4 pulses at different frequencies, possibly include models with 6 s of constant stimulation at different CamCa4 concentrations

Next Steps:
4. Update fitness function for experimental data
5. Specify parameters to vary in a separate file, 
     possibly with xml specifications similar to the model instead of xpath - see neurord3_params_sims.py and neurord_2params_sims.py

6a. add features to fitness function, e.g. peak amplitude and width - need to specify stimulus onset time, as with injection current 
6b. add fitness functions to file in ajustador.  Attach fitness values to fit object

7. Add variation of IC parameters
'''

import ajustador as aju
import numpy as np
import loadconc
import tempfile
from ajustador import drawing
import glob

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='camkii/'  #location of data and model files.
#Set of model files that have first part of file name in common.  All included files must be in same directory.
model_set=dirname+'Model-CKnew-Cahz'
exp_set=dirname+'Model-CKold-Cahz' #set of data files corresponding to model files; files may contain several molecules
mol=['CKpCamCa4'] #which molecule(s) to match in optimization
tmpdir='/tmp/out' 

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=10
popsize=8

P = aju.xml.XMLParam
#list of parameters to change/optimize
params = aju.optimize.ParamSet(P('CK2_fwd_rate', 0, min=0, max=1e-6, xpath='//Reaction[@id="CKCam_pow2"]/forwardRate'),
                               P('CK4_fwd_rate', 0, min=0, max=1e-9, xpath='//Reaction[@id="CKCam_pow4"]/forwardRate'),
                               P('CK3_fwd_rate', 0, min=0, max=1e-12, xpath='//Reaction[@id="CKCam_pow3"]/forwardRate'),
                               P('CK1_CKp2_fwd_rate', 0, min=0, max=1e-9, xpath='//Reaction[@id="CK1_CKpCam_pow2"]/forwardRate'),
                               P('CK2_CKp1_fwd_rate', 0, min=0, max=1e-9, xpath='//Reaction[@id="CK2_CKpCam_pow1"]/forwardRate'),
                               P('CK2_CKp2_fwd_rate', 0, min=0, max=1e-12, xpath='//Reaction[@id="CK2_CKpCam_pow2"]/forwardRate'))

###################### END CUSTOMIZATION #######################################
exp = aju.xml.NeurordResult(exp_set)

###################### fitness_functions - move to nrd_fitness.py in ajustador #######################################
def nrd_output_conc(sim_output,species):
    #may need to add specification of trial and/or voxel
    pop1count = sim_output.counts().xs(species,level=2)
    volumes,PUVC=sim_output.volumes()
    tot_vol=np.sum(volumes)
    pop1conc=pop1count.sum(axis=0,level=1)/tot_vol/PUVC  #sum across voxels, level=0 sums across time
    return pop1conc

def specie_data_concentration_fitness(*, voxel=0, species_list, trial=0):
    def fitness(sim, measurement, full=False):
        #2 is level of multi-index for species, 0 is level of multi-index for voxel
        fitarray=np.zeros(len(species_list))
        diffarray={}
        for i,species in enumerate(species_list):
            pop1conc=nrd_output_conc(sim.output[0],species)
            wave1y=pop1conc.values[:,0]
            wave1x=pop1conc.index
            pop2 = measurement.waves[species].wave
            # Note: np.interp(x1,x2,y2) returns values for y2 corresponding to x1 timepoints
            pop1y=np.interp(pop2.x,wave1x,wave1y)
            diff = pop2.y - pop1y
            fitarray[i]=float((diff**2).mean()**0.5)
            diffarray[species]=diff
        #print('fit',species_list,fitarray)
        fitness=np.mean(fitarray)
        if full:
            return diffarray
        else:
            return fitness
    return fitness

def specie_sim_concentration_fitness(*, voxel=0, species_list, trial=0):
    def fitness(sim, measurement, full=False):
        print('fitness', sim,measurement)
        fitarray=np.zeros((len(species_list),len(sim.output)))
        fit_dict={}
        for i,species in enumerate(species_list):
            fit_dict[species]={}
            for j,stim_set in enumerate(sim.output):
                print ('aligned?', stim_set.file.filename, stim_set.injection, measurement.output[j].file.filename, measurement.output[j].injection)
                pop1=nrd_output_conc(stim_set,species)
                #pop1 = stim_set.concentrations().loc[voxel, :, species, trial]
                if isinstance(exp,aju.xml.NeurordResult):
                    #pop2 = measurement.output[i].concentrations().loc[voxel, :, species, trial]
                    pop2 = nrd_output_conc(measurement.output[j],species)
                #else - do stuff with waves
                diff = pop2 - pop1
                fit_dict[species][stim_set.injection]=float((diff**2).mean()**0.5)
                fitarray[i][j]=float((diff**2).mean()**0.5)
        fitness=np.mean(fitarray)
        print (fitarray)
        if full:
            return fit_dict
        else:
            return fitness
    return fitness

fitness = specie_sim_concentration_fitness(species_list=mol)

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

########################################### Done with fitting

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
print(fit.param_names())
print(fit.params.unscale(fit.optimizer.result()[0]))
print(fit.params.unscale(fit.optimizer.result()[6]))

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

#to plot data (clicking on point in fit history doesn't work)
#if other than last sim is desired, assign that value to fitnum
def plot_traces(fitX,mollist,fitnum=-1):
    from matplotlib import pyplot as plt
    import os
    plt.ion()
    bestdir=fitX[fitnum].tmpdir.name
    stim_set=glob.glob(bestdir+'/model*.h5')
    print ('PLOTTING', stim_set, 'parameters: ',fit[fitnum])
    fig,axes=plt.subplots(len(mollist),1,sharex=True,figsize=(6,9))
    fig.canvas.set_window_title(fitX.model)
    axis=fig.axes
    for j,fname in enumerate(stim_set):
        best=aju.nrd_output.Output(fname)
        label=os.path.basename(aju.nrd_output.Output(fname).file.filename)[0:-3]
        for i,mol in enumerate(mollist):
            simdata=nrd_output_conc(best,mol)
            axis[i].plot(simdata.index,simdata.values[:,0],label=label)
            #if isinstance(exp,aju.xml.NeurordResult):
            expdata=nrd_output_conc(exp.output[j],mol)
            axis[i].plot(expdata.index,expdata.values[:,0])
            #else
            #axis[i].plot(exp.waves[mol].wave.x,exp.waves[mol].wave.y)
            axis[i].set_ylabel(mol)
        axis[i].set_xlabel('time, msec')
    axis[i].legend()

fitnum=-1
plot_traces(fit,mol)

