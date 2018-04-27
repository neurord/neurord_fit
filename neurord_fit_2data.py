#How to use :doc:`ajustador` to fit a NeuroRD model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fits a single reaction model (Model_pSynGap) to data.
'''

import ajustador as aju
import numpy as np
import loadconc
import tempfile
from ajustador import drawing
import glob    

#model is the xml file that contains the neuroRD model to simulate and adjust parameters
dirname='.'  #where is data stored.  Multiple datafiles allowed
model='Model_syngap_ras.xml'
exp_name='walkup_JBC_2015' #name of data file selected from dirname; each file may contain several molecules
mol=['pSynGap','RasGDP'] #which molecule(s) to match in optimization
tmpdir='/tmp/'+exp_name

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

#exp is the data to fit  
csvs = sorted(glob.glob('{}/*.csv'.format(dirname)))
allexp = [loadconc.CSV_conc(name) for name in csvs]
allexp={series.name.split('/')[-1].split('.')[0]:series for series in allexp}
exp=allexp[exp_name]

###################### fitness_functions - move to nrd_fitness.py in ajustador #######################################
def nrd_output_conc(sim_output,species):
    #may need to add specification of trial and/or voxel
    pop1count = sim_output.counts().xs(species,level=2)
    volumes,PUVC=sim_output.volumes()
    tot_vol=np.sum(volumes)
    pop1conc=pop1count.sum(axis=0,level=1)/tot_vol/PUVC  #sum across voxels, level=0 sums across time
    return pop1conc

def specie_concentration_fitness(*, voxel=0, species_list, trial=0):
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

fitness = specie_concentration_fitness(species_list=mol)

############ Test fitness function
#sim = aju.xml.NeurordSimulation('/tmp', model=model, params=params)
#cp /tmp/???/model.h5 modelname.split('.')[0]+'.h5'
#sim2=aju.xml.NeurordResult('Model_syngap_ras.h5')
#print(fitness(sim2, exp))
################

fit = aju.optimize.Fit(tmpdir, exp, model, None, fitness, params,
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

'''
NeurordSimulation(<TemporaryDirectory '/tmp/walkup_JBC_2015/tmp6n3kec2r'>, phos_fwd_rate=Param phos_fwd_rate=0.9670002607916861 phos_rev_rate=Param phos_rev_rate=0.21662546285518958 phos_kcat_rate=Param phos_kcat_rate=0.014834835304915908 gap_kf_rate=Param gap_kf_rate=0.33152534037329834 gap_kb_rate=Param gap_kb_rate=0.10923292188031208 gap_kcat_rate=Param gap_kcat_rate=0.5953407411944236 Pgap_kf_rate=Param Pgap_kf_rate=0.11143400862714481 Pgap_kb_rate=Param Pgap_kb_rate=0.9784048303259848 Pgap_kcat_rate=Param Pgap_kcat_rate=0.828326863711133)

vs
 
>>> fit.param_names()
['phos_fwd_rate', 'phos_rev_rate', 'phos_kcat_rate', 'gap_kf_rate', 'gap_kb_rate', 'gap_kcat_rate', 'Pgap_kf_rate', 'Pgap_kb_rate', 'Pgap_kcat_rate']
>>> print(fit.params.unscale(fit.optimizer.result()[0]))
[0.00038412011773427307, 0.22458737290234021, 0.2204795879164167, 0.24461004409872056, 0.2414575433038656, 0.095997833187105863, 0.12906754484253635, 0.092105300850746502, 0.26312214582883092]
'''

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

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

