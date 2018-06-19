How to use :doc:`ajustador` to fit a NeuroRD model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. clone the following packages:

   a. neurord/ajustador
   b. neurord/neurord_fit
      
2. get neurord-3.2.3-all-deps.jar from neurord/stochdiff/releases
3. Get into python3 as follows:
   
   PYTHONPATH=$PYTHONPATH:/full/path/to/ajustador/:/full/path/to/neurord_fit/ python3
   
4. Run the example commands shown in neurord_fit_simple.py. 
   
A. To use this to optimize rate constants for other reaction files, do the following steps:

1. create NeuroRD model with the reaction (or reactions) to optimize
2. update the parameters to optimize (params= ...)
3. To optimize using real data (or a different model), csv formatted file with molecule concentration or percent change.
   
   a. data format
      
      - 1st column is time, either ms, sec or min.  If sec or min, indicate units separated by space from "Time"
	
      - additional columns contain concentration or FRET.  header is name of molecule to compare to in simulation
	
	+ If concentration is mMolar or uMolar, indicate units separated by space from molecule name (default is nanoMolar)
	  
	+ if Fret, can scale concentration by a user specified value.  Specify units as % and indicate scale factor separated by space from units
	  
   b. specify
      
          - norm_method='percent'
	  
	  - start_stim=<value>
	  
      for matching simulations to percent change, and indicating the time of stimulation.
      The time between 0 and stim_start will be used for baseline (denominator for percent change)
      
   c. use exp=loadconc.CSV_conc_set(exp_name,stim_time=start_stim) for loading experimental data, default stim_time is 0  
      
   d. specify fitness=nrd_fitness.specie_concentration_fitness(species_list=mol,start=start_stim,norm=norm_method).
      - in fitness function, default normalization divides the difference between model and data by the peak
      
4. fitness function can match a set of molecules - specify mol=['mol1','mol2']
5. Optimization can match a set of different model and experimental files.  Set of files should have 1st part of filename in common, and differ only by the suffix.  The model files need to have matching suffixes, e.g.
   
   - model_stim1.xml, model_stim2.xml, exp_stim1.csv, exp_stim2.csv

B. Examples in the repo

1. neurord_fit_simple.py: fits a single simulation with one reaction to a previous simulation, matches one molecule concentration
2. neurord_fit_2data.py: fits a single simulation to a single data file, matches two different molecule concentrations
3. neurord_fit_CamKII.py: fits a set of simulations to a previous set of simulations
4. neurord_fit_CamKIIdata.py: fits a set of simulations to data on CamKII autophorphorylation
5. neurord_fit_CamKII_PP1.py: fits dephosphorylation parameters to achieve a 1 uM basal phopshoCamKII
6. neurord_fit_fret.py: fits a set of simulations to a set of experimental data, normalizes simulation by baseline (norm_method='percent', start_stim=100) to match percent change in FRET in data

C. Planned improvements to optimization

1. possibly specify parameters to vary using xml specs similar to the model instead of xpath 
2. Allow variations to kcat while constraining Km for enzyme reactions, or vary kb while constraining KD for bimolecular rxn, or total concentration
3. Use features in the fitness function, such as peak amplitude, width, decay time constant 
4. allow stimulation onset of model and experimental data to differ, and align the waveforms for the fitness function
