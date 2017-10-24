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
2. create an array of value versus time with format the same as exp.output.counts().loc[0, :, 'glu', 0]
3. update pop2= in the fitness definition
4. update the pop1 - pop2 comparison to match time samples
