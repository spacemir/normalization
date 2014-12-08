To run simulations NEST 2.2.2, pyplot, numpy and matplotlib needs to be installed.


BASIC STRUCTURE
HC_ModelA is base class describing the structure of the hypercolumn model
HC_STD inherits from HC_ModelA and define new build and connect methods where short-term
depression is added

hc_utils contains methods that simulate hypercolumn activity and calculates e. g. statistics

hc_visualization, simple_plots and utility_functions are all needed for calculating and
visualizing test results.  

"run_simulation" is a script by which a certain simulation is ran, first command line arguments
give which test, second if the test is ran from a previously saved file.

"h"= hypercolumn test. The network activity for a single input vector is studied,
	by plotting spike histograms, raster plots, voltage traces, and average spike frequencies for 
	neurons in each minicolumn and the basket cell population
"io" = IO test
"fir" = FIR test
"2mc" = Test where the output from a network with two active hypercolumn is plotted as 
	both average output, and output relations, for different inputs to the two minicolumns

To run e. g. the io test thus give the following terminal command:
> python run_simulation io 


The different files parameters_xxx contain parameter dictionaries, used in the run_simulation script,
to use another parameter dictionary this has to be changed in this file.

The same is true for switching between the standard hypercolumn class "HC" and "HC_STD"




