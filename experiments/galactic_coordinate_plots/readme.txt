OVERVIEW

The experiment in this folder is to create the same plots as Ibata 2020.
Ibata 2020 created color plots with galatic latitude and longitude. 
The colorplots parameters of the color plots are heliocentric distance, line of sight velocity, temperature
We will do heliocentric distance and angular momentum 


INPUTS 
(I0) - ../simulations/outputDATA/*galacticmodel*/*globular_potential*/steams/
(I1) - mini_list_of_clusters.txt
(I2) = ../../functions/clean_streams

OUTPUTS
a series of png images named
(O0) coordinateSystem_galacticPotential_GCPotential_**NAME_OF_QUANTITY**.png
	EX. galactic_coordinates_PII_Plummer_D_NGC5024.png
	Ex. galactic_coordinates_PII_Plummer_recent_1_Gyrs_D_NGC5986.png
	galactic_coordinates_PII_Plummer_recent_1_Gyrs_PMB_short_color_bar__NGC6093.png
(O1) coordinateSystem_galacticPotential_GCPotential_**NAME_OF_QUANTITY**_nClusters_NCLUSTERS.png
	ex. galactic_coordinates_PII_Plummer_D_nClusters_7
(O2) output_heliocentric_profile/*
	Portions of profiles saved into its own directory
(O3) output_heliocentric_profile_5_kpc/*
(O4) output_clean_streams_helio_3_kpc/* 
(O5) output_heliocentric_images_gifs/*
	Gifs created out of profiles using online free software
(O6) density_plots/
	
CODE
(C0) - galactic_latitude_longitude_color_plot
(C1) - galactic_coordinate_color_plot
(C2) - galactic_plot_individual.py
(C3) - galactic_plot_multiple.py
(C4) - galactic_profiles.py
(C5) - density_plots.py
(C6) - concatenate_profiles.py

the code file (C0) grabs data from (I0). This code creates colorplots in the galactic latitude and longitude. The colors are quantities such as the heliocentric distance, and galactic z angular momentum. The user must specify which galatic model to use, for example PII< in the code. The code will then search in the /galacticmodel/steams/ folder and grab each h5 data. It will plot the desired quantity for each star in the steams folder from each h5 file. 

(C1) does the same but in galactic centric coordaintes 

(C2) A combination of C0 and C1. I.e. the user can just change arguments to create a plot in either equatorial, galactic, or galactocentric coordaintes. This has outputs of form of 

(C3) Plots multiple streams onto the same plot. The user can choose the coordinate system and the color bar. Also, (I1) is the argument for chooseing which clusters to plot. The user must include the cluster name and the number of years to concatenate. 5 will include all, where 1 would only include the last one billion years 

(C4) Plots profiles of streams. I.e. only windows based on helio-centric dstance. 
The user must edit the arguments down_lim_list and up_lim_list to change the width of each window and the step size between windows. 

(C5) Creates plots of the density of each GCS ejected stellar material. 
We plot this is either units of counts per KPC or counts per stradian. The outputs are in (O6) density_plots. 
I personally think this is one of the best written scripts I have

(C6) makes the trifecta images in density plots.


IF TAKING THESE FROM GITHUB, SOME OF THESE OUTPUTS MIGHT HAVE BEEN SURPRESSED
