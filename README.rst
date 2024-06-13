How to use :doc:`ajustador` to fit a NeuroRD model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. clone the following packages:

   a. neurord/ajustador
   b. neurord/neurord_fit
      
2. get neurord-3.3.0-all-deps.jar from neurord/stochdiff/releases
3. Get into python3 as follows:
   
   PYTHONPATH=$PYTHONPATH:/full/path/to/ajustador/:/full/path/to/neurord_fit/ python
   
4. Run the example commands shown in neurord_fit_simple.py. 
   
   A. To use this to optimize rate constants for other reaction files, do the following steps:

      1. create NeuroRD model, which includes reactions and concentrations to optimize, as well as morphology and (optionally) stimulation
         - For computational efficiency, having one or very few mesh elements is recommended
         - see https://github.com/neurord/stochdiff#readme for more information on model specification

      2. update the parameters to optimize (params= ...)

         a. without any constraints

            - list the molecule name, which is user defined, starting value, maximum value, minimum value, and the xpath in that order

               + max and min values are optional

            - example:

                  P = aju.xml.XMLParam
                  params = aju.optimize.ParamSet(P('Racact_fwd_rate',2.79888e-07, min=2.8e-9, max=2.8e-5, xpath='//Reaction[@id="RacGDP+Kal--pKalRacGDP"]/forwardRate'))

         b. with multiplicative constraints

            - if two parameters relate to each other (i.e. rate A is X times rate B), then make X an integer stored in the variable constant

            - make rate B a string stored in the variable fixed

            - list the molecule name, starting value, fixed, constant, and the xpath in that order

            - example:

               P = aju.xml.XMLParam
               params = aju.optimize.ParamSet(P('Racact_bckd_rate',0.0001658, fixed='Racact_kcat_rate',constant=4, xpath='//Reaction[@id="RacGDP+Kal--pKalRacGDP"]/reserveRate'))

         c. with summative constraints

            - if you have parameters that sum together (i.e. molecule A equals X minus molecule B), then make X an integer stored in the variable constant

            - make a dictionary that has two entries, one being molecules and one being radius

               + the entry molecules should have a list of the molecules being summed together given as strings

                  - if the molecule is a surface density, it must have an ending of "_dens"

               + the entry radius should have the radius of the cell being modeled given as an integer
               
            - list the molecule name, starting value, fixed, constant, and the xpath in that order

            - example:

               P = aju.xml.XMLParam
               params = aju.optimize.ParamSet(P('CamCa4', 1000, fixed={'molecules':['AC1CamATP_dens'],'radius':1.5}, constant=1000, xpath='//NanoMolarity[@specieID="CamCa4"]'))
               
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
      
      4. fitness function can match a set of molecules, each of which can be the sum of molecule species
		- specify the molecules as a dictionary:mol={'mol1':['mol1'],'mol2':['subspeciesA', 'subspeciesB']}
		- the key matches the data to optimize (e.g. the name of the molecule in the csv data file)
		- the value is a list of subspecies to sum

      5. Optimization can match a set of different model and experimental files.  Set of files should have 1st part of filename in common, and differ only by the suffix.  The model files need to have matching suffixes, e.g.
   
         - model_stim1.xml, model_stim2.xml, exp_stim1.csv, exp_stim2.csv
   
      6. After running the optimization, the parameter values are in the file outcmaesrecentbest.dat, which also contains the fitness value. The order of parameters is the same as specified in the params= ... line in the neurord_fit.py file. The best parameter values are in the last line in outcmaesrecentbest.dat, or the line with the smallest fitness value.

   B. Examples in the repo

      1. neurord_fit_simple.py: fits a single simulation with one reaction to a previous simulation, matches one molecule concentration
      2. neurord_fit_2data.py: fits a single simulation to a single data file, matches two different molecule concentrations
      3. neurord_fit_CamKII.py: fits a set of simulations to a previous set of simulations
      4. neurord_fit_CamKIIdata.py: fits a set of simulations to data on CamKII autophorphorylation
      5. neurord_fit_CamKII_PP1.py: fits dephosphorylation parameters to achieve a 1 uM basal phopshoCamKII
      6. neurord_fit_fret.py: fits a set of simulations to a set of experimental data, normalizes simulation by baseline (norm_method='percent', start_stim=100) to match percent change in FRET in data

   C. Planned improvements to optimization

      1. Use features in the fitness function, such as peak amplitude, width, decay time constant 
      2. allow stimulation onset of model and experimental data to differ, and align the waveforms for the fitness function
