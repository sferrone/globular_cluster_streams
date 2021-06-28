OVERVIEW
This experiment comes from analyzing the colorbars of each stellar stream in Ibata et al 2021 and creating a list of constraints for each stream. We then check if the orbits of the globular clusters pass through this interval. If it does, we save it and when it passes in json files. We then develop other code to plot this phenomena



CODE/
(C0)	- create_interval_json.py
(C2)	- galactic_coordiantes_with_stream.py
(C3)	- generate_all_json_initial_coniditions.py
(C4)	- generate_subset_json_initial_conditions.py
(C5)	- intersection_json.py
(C6)	- plot_all_matches.py
(C7)	- plot_match_stream_and_cluster.py
(C8)	- run_local_matching.sh
(C9)	

INPUTS 
(I0)	- mini_list_of_clusters.txt
(I1) the_orbits/
	- a folder containing .dat files of the orbits of the globular clusters forward in 
	- This will probably be in a different place depending on the machine

OUTPUTS 
(O1) OUTPUTS/

(O1) output_stream_to_cluster/



EXPLAINATIONS 

(C0) is a function used to create the JSON objects containing the initial conditions for each stream.

(C3) essentially a tabulated list of the constrains for each stellar stream. This code executes C0 to initialize all json objects and stores them in outputs/. This code was only updated as a mistake was found, but it was never compiled after its first execution.

(C4) We only update the initial conditions of a couple streams. This way we don't rewrite all the json objects. 

(C5) This code is intense. It checks the orbital evolution of each globular cluster to see which satisfy the constrains of each stream specified by the JSON objects. If yes, the JSON object is appended with the globular cluster. This code will check each of the streams according to I0. It uses the .dat files in (I1). 

(C6)	This code plots the orbital evolution of each globular cluster who satisfy the criteria of each stellar stream. The JSON objects of the (O0) outputs/ are checked. These outputs are saved in a folder system under (O1). 

(C7)  	Honestly this code is trash. It tries to plot all the potential GCs of the streams onto one plot. It is saved into stream_matches/

(C2)	This is to do direct checks between a globular cluster and a stream 

(C8) executes the matching 