# Code used in Pisokas, Heinze, Webb 2019, Biorxiv and Pisokas, Heinze, Webb (2020), eLife. 

This is the code accompanying the Pisokas, Heinze, Webb (2019), Biorxiv and Pisokas, Heinze, Webb (2020), eLife. It uses the code in the directory Kakaria_deBivort_2017/ for running the simulations. The code in Kakaria_deBivort_2017/ is modified version of the code made available by Kakaria and de Bivort, 2017.

The code in Exp_Recurrent_PB_EB/ is the code developed for the paper Pisokas, Heinze, Webb (2019), Biorxiv and Pisokas, Heinze, Webb (2020), eLife.


# Usage

Both directories and subdirectories must be added in Matlab path. 

Running the script plot_explore_bump_jump_with_P_EGs.mlx in Matlab2016b will reproduce the plots included in the paper (the code has been tested with Matlab2016b and Matlab2020a). 

Running the script make_plots.Rmd using R will reproduce the plots of response time and heading change rate plots in the paper. 

Running the script collect_stats_long_run.m will run the simulations and collect all the data used in the plots (this code has been tested with Matlab2016b). The script contains instructions and list of the ways it was invoked to procude each of the data sets. 


# References
Kakaria KS, de Bivort BL. Ring Attractor Dynamics Emerge from a Spiking Model of the Entire Protocerebral Bridge. Frontiers in Behavioral Neuroscience. 2017; 11:8. doi: 10.3389/fnbeh.2017.00008.
