#How to use :doc:`ajustador` to fit a NeuroRD model to FRET data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''This demonstration fiting a model of signaling pathways from norepinephrine to Epac_FRET and PKA
 to FRET data with two difference conditions, in this case diameter of structure.
 This also shows optimizing both initial conditions and reaction reates
Note that stimulation must occur at same time in simulation and experiment
Necessary to specify time of stimulation, in seconds, both in this file and in Model stimulation file
xEliminate calcium & calmodulin reactions to speed simulation, and decrease num comps in 6.5u from 9 to 5.
Need to tweak the conversion from model Epac1cAMP to data FRET signal else fitness function won't work.
'''

import ajustador as aju
import numpy as np
from ajustador import drawing,loadconc,nrd_fitness
from ajustador.helpers import converge,save_params
import os

dirname='fret_cAMP/'  #where data and model file are stored.  Can be different than current directory. Multiple datafiles allowed
#model_set is the set of xml files that contains the neuroRD model to simulate and adjust parameters
#Note that larger volumes take longer to simulate, so use small depth2D for optimization
model_set='Model_VinceCaPKAsubset_ACsame-dia'
exp_name='FretPercent' #name of data file selected from dirname; each file may contain several molecules
mol=['Epac1cAMP'] #which molecule(s) to match in optimization
tmpdir='/tmp/AC_PDEsame_'+dirname
start_stim=140  #time of onset of stimulation, in seconds; should be able to extract from model files; stim time must match data
norm_method='percent' #convert molecule concentration into a percent change from baseline.  Use for comparing to FRET data

# number of iterations, use 1 for testing
# default popsize=8, use 3 for testing
iterations=50
popsize=6 # reduce from 8 to avoid saturating processors for longer neurord simulations
test_size=25

os.chdir(dirname)
#this command indicates that experimental data are concentration in csv formatted files
exp=loadconc.CSV_conc_set(exp_name)

P = aju.xml.XMLParam
#list of parameters to change/optimize
#improvement: allow rate (multiply all kf,kb,kcat) or KM (vary what? Kb?)
params = aju.optimize.ParamSet(P('phospLRGs_fwd_rate', 0.00125e-3, min=0.0001e-3, max=0.1e-3, xpath='//Reaction[@id="phospLRGs1"]/forwardRate'),
                               P('phospLRGs_bak_rate', 0.16e-3, min=0.01e-3,max=1e-3,xpath='//Reaction[@id="phospLRGs1"]/reverseRate'),
                               P('phospLRGs_cat_rate', 0.04e-3, min=0.001e-3,max=1e-3,xpath='//Reaction[@id="phospLRGs2"]/forwardRate'),
                               P('dphospD1R_fwd_rate', 0.003e-3, min=0.0001e-3,max=0.1e-3,xpath='//Reaction[@id="dphospD1R"]/forwardRate'),
                               P('GasGTP_hydrolysis', 0.5e-3, min=0.01e-3,max=5e-3,xpath='//Reaction[@id="GasGTP_disso"]/forwardRate'),
                               P('AC1_GasGTP_GAP', 1e-3, min=0.01e-3,max=10e-3,xpath='//Reaction[@id="AC1_GasGTP_GAP"]/forwardRate'),
                               P('dphospPDE4_fwd_rate', 0.00625e-3, min=0.0001e-3,max=0.1e-3,xpath='//Reaction[@id="dphospPDE4"]/forwardRate'),
                               P('AC1_conc', 1000, min=1,max=10e3,xpath='//PicoSD[@specieID="AC1"]'),
                               P('AC1CamATP_conc', 200, min=1,max=1000,xpath='//PicoSD[@specieID="AC1CamCa4ATP"]'),
                               P('AC5_conc', 200, min=1,max=1000,xpath='//PicoSD[@specieID="AC5"]'),
                               P('PDE4_conc', 1000, min=1,max=10e3,xpath='//NanoMolarity[@specieID="PDE4"]'),
                               P('pPDE4_conc', 50, min=1,max=1000,xpath='//NanoMolarity[@specieID="pPDE4"]'))

###################### END CUSTOMIZATION #######################################
fitness = nrd_fitness.specie_concentration_fitness(species_list=mol,start=start_stim,norm=norm_method)
fit = aju.optimize.Fit(tmpdir, exp, model_set, None, fitness, params,
                       _make_simulation=aju.xml.NeurordSimulation.make,
                       _result_constructor=aju.xml.NeurordResult)
fit.load()
#fit.do_fit(iterations, popsize=popsize,sigma=0.3)
fit.do_fit(iterations, popsize=popsize,seed=62839)
mean_dict,std_dict,CV=converge.iterate_fit(fit,test_size,popsize)

########################################### Done with fitting, look at results

#to look at fit history
aju.drawing.plot_history(fit,fit.measurement)

#to look at centroid [0] or stdev [6] of cloud of good results:
#to recall names of parameters
for i,p in enumerate(fit.params.unscale(fit.optimizer.result()[0])):
    print(fit.param_names()[i],'=',p, '+/-', fit.params.unscale(fit.optimizer.result()[6])[i])

save_params.save_params(fit,0,1)
''' Results
1. Convergence reached, though variability relatively large, and decay phase not good. 
Perhaps need to adjust convergence criteria to normalize to mean fitness.
                                                                       previous values:
phospLRGs_fwd_rate = 1.49376180957e-06 +/- 9.25145585683e-08           0.00125e-03 - similar
phospLRGs_bak_rate = 7.78327479984e-05 +/- 4.95733956741e-06           0.16e-03  - opt 2x smaller 
phospLRGs_cat_rate = 5.30217929784e-05 +/- 7.7488161955e-07            0.04e-03 - opt similar
dphospD1R_fwd_rate = 1.12592920188e-07 +/- 6.85755643513e-08           0.003e-03 - opt 27x smaller
GasGTP_hydrolysis = 0.000460175493085 +/- 4.70658013805e-06            0.5e-3 - opt similar
AC1_GasGTP_GAP = 0.000537311562671 +/- 8.80953846399e-05               1e-3 - opt is half value
dphospPDE4_fwd_rate = 7.58695147211e-06 +/- 4.88583688336e-08          0.00625e-03 - opt is 10x smaller
smaller phospLRGs dephos rate will produce net greater R phosphorylation - verify by evaluating files in /tmp/fret
accounted for by lower GasGTP GAP activity?
10x lower dephos PDE4 - should see great pPDE4 - verify by evaluating files in /tmp/fret

2. Use this test set for optimization of initial conditions.  Allow both AC and PDE4 to vary.
phospLRGs_fwd_rate = 1.35727020376e-06 +/- 1.80794893966e-07
phospLRGs_bak_rate = 0.000153244527315 +/- 2.15329689388e-05
phospLRGs_cat_rate = 4.48123423014e-05 +/- 2.1570514001e-06
dphospD1R_fwd_rate = 4.17210061666e-07 +/- 2.67090322749e-07
GasGTP_hydrolysis = 0.000194143562613 +/- 2.24639192471e-05
AC1_GasGTP_GAP = 0.000407303048865 +/- 0.000123488151497
dphospPDE4_fwd_rate = 3.6303256312e-06 +/- 2.18413949353e-07
AC1_conc = 3938.39760923 +/- 250.54612832
AC1CamATP_conc = 159.916164109 +/- 21.1579947125
AC5_conc = 504.867439242 +/- 23.3965158586
PDE4_conc = 4347.9687042 +/- 258.534274792
pPDE4_conc = 56.4948995096 +/- 1.98002521843

Good fit for 1.5u, but not  for 6.5u - this would be an AC same + PDE4 same simulation
>>2018jun4_fretpercent_optAC_PDE.png
seed2 (june7)
phospLRGs_fwd_rate = 3.97074601743e-06 +/- 6.49225175073e-07
phospLRGs_bak_rate = 0.000502994563306 +/- 7.72799593502e-05
phospLRGs_cat_rate = 2.48348696394e-05 +/- 8.78336451609e-06
dphospD1R_fwd_rate = 4.78308505489e-06 +/- 8.23862644244e-07
GasGTP_hydrolysis = 0.000482153787693 +/- 7.56873501461e-05
AC1_GasGTP_GAP = 0.000446300058826 +/- 0.00048969685802
dphospPDE4_fwd_rate = 3.47524283979e-07 +/- 8.22241832131e-07
AC1_conc = 1380.34167736 +/- 579.314220409
AC1CamATP_conc = 809.523109891 +/- 108.041041992
AC5_conc = 332.577281536 +/- 67.8872680684
PDE4_conc = 9050.94955877 +/- 905.915261881
pPDE4_conc = 113.235558268 +/- 7.10899862915

3. To constrain AC differently for 1.5 and 6.5, then use PDE4same, and only tune PDE4
>>>>>>Not as good.  saved images
phospLRGs_fwd_rate = 1.2035426655e-06 +/- 7.14283853283e-08
phospLRGs_bak_rate = 0.000321277059461 +/- 9.52279769966e-06
phospLRGs_cat_rate = 4.00856246999e-05 +/- 7.76624500839e-07
dphospD1R_fwd_rate = 1.47244215222e-06 +/- 1.0531311591e-07
GasGTP_hydrolysis = 0.000502451979063 +/- 7.70030047769e-06
AC1_GasGTP_GAP = 0.00084557569035 +/- 7.63560640647e-05
dphospPDE4_fwd_rate = 7.02464993607e-06 +/- 1.03653799369e-07
PDE4cAMP_conc = 283.201328172 +/- 7.17804102159
PDE4_conc = 1175.35399912 +/- 61.4213870296
pPDE4_conc = 54.6653064118 +/- 0.948036317897

ACdiff->seed2 (june6)
phospLRGs_fwd_rate = 7.9978158674e-07 +/- 1.74754651971e-07
phospLRGs_bak_rate = 0.00012603976567 +/- 1.5880911789e-05
phospLRGs_cat_rate = 3.00561717756e-05 +/- 1.62507149005e-06
dphospD1R_fwd_rate = 2.38040186357e-06 +/- 1.69974533172e-07
GasGTP_hydrolysis = 0.000734972557666 +/- 1.34425813135e-05
AC1_GasGTP_GAP = 0.000598949149457 +/- 0.00010844230616
dphospPDE4_fwd_rate = 5.38076632634e-06 +/- 1.60873390266e-07
PDE4cAMP_conc = 63.2549614948 +/- 17.9041015299
PDE4_conc = 1604.74634895 +/- 142.238802519
pPDE4_conc = 70.8118670051 +/- 1.71918861

Do results support the AC/PDE4 ratios used for hand-tuning?

3a. Add ability to tune kcat without changing KM - will need to identify the correct binding reaction and modify those
3b. Add ability to constrain Kb >= 4kcat
Use this model to test above.
'''
