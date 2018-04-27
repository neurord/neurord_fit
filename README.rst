How to use :doc:`ajustador` to fit a NeuroRD model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. clone the following packages:

   a. neurord/ajustador
   b. neurord/neurord_fit
      
2. get neurord-3.2.3-all-deps.jar from neurord/stochdiff/releases
3. Get into python3 as follows:
   PYTHONPATH=$PYTHONPATH:/full/path/to/ajustador/:/full/path/to/neurord_fit/ python3
4. Run the example commands shown in neurord_fit.py.  Make sure to update the path to the model (model = aju.xml.open_model('neurord_fit/Model_simple.xml') and output data as appropriate.
   
To use this to optimize rate constants for other reaction files, do the following steps:

1. create Model with the reaction (or reactions) to optimize
2. update the parameters to optimize
3. To optimize using real data (or a different model), create an array of value versus time with format the same as exp.output.counts().loc[0, :, 'glu', 0].  I.e., write a new function to read in text or csv files.
4. interpolate pop1 or pop2 so that the values of each correspond to the same samples in time
5. Optionally create an entirely different fitness function
6. Expand fitness function and simulation control to run multiple timulation files / conditions, and evaluate the result of all of them to determine the fitness
