How to use :doc:`ajustador` to fit a NeuroRD model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. clone the following packages:

   a. neurord/ajustador
   b. neurord/neurord_fit
      
2. get neurord-3.2.3-all-deps.jar from neurord/stochdiff/releases
3. Get into python3 as follows:
   PYTHONPATH=$PYTHONPATH:/full/path/to/ajustador/:/full/path/to/neurord_fit/ python3
4. Run the example commands shown in neurord_fit.py. 
   
To use this to optimize rate constants for other reaction files, do the following steps:

1. create NeuroRD model with the reaction (or reactions) to optimize
2. update the parameters to optimize (params= ...)
3. To optimize using real data (or a different model), csv formatted file with molecule concentration or percent change
   a. use exp=loadconc.CSV_conc_set(exp_name) for loading experimental data
   b. specify
          norm_method='percent'
	  start_stim=<value>
      for matching simulations to percent change, and indicating the time of stimulation.
      The time between 0 and stim_start will be used for baseline (denominator for percent change)
   c. specify fitness=neurord_fit.specie_concentration_fitness(species_list=mol,start=start_stim,norm=norm_method)
4. fitness function can match a set of molecules - specify mol=['mol1','mol2']
5. Optimization can match a set of different model and experimental files.  Set of files should have 1st part of filename in common, and differ only by the suffix.  The model files need to have matching suffixes, e.g.
   model_stim1.xml, model_stim2.xml, exp_stim1.csv, exp_stim2.csv

Examples in the repo
1. neurord_fit_2data.py: fits a single simulation to a single data file, matches two different molecule concentrations
2. neurord_fit_CamKII.py: fits a set of simulations to a previous set of simulations
3. neurord_fit_fret.py: fits a set of simulations to a set of experimental data, normalizes simulation by baseline (norm_method='percent', start_stim=100) to match percent change in FRET in data

Planned improvements to optimization
1. Allow optimization of initial conditions
2. possibly specify parameters to vary using xml specs similar to the model instead of xpath (see neurord3_params_sims.py and neurord_2params_sims.py for examples)
3. add features to the fitness function, such as peak amplitude, width, decay time constant (Attach fitness values to fit object?)
4. specify stimulation onset of model and experimental data separately, and align the waveforms for the fitness function
