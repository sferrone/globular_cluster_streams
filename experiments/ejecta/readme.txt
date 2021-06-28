Overview
This experiment is supposed to show information about when and where the stars escape

Code files
(C1) code/animate_fugitives_rz.py
(C2) code/animate_fugitives_xy.py
(C3) code/escaper_histograms.py
(C4) code/print_clusters.sh

Inputs
(I1) list_of_clusters.txt

Outputs
(O1) outputs/GCNAME_galacticPotential_GCPotential_fugitives_rz.gif 
    Example:    NGC6528_PII_Plummer_fugitives_rz.gif
(O2) outputs/escaper_histogram_galacticPotential_GCPotential_GCNAME.png 
    Example: escaper_histogram_PII_Plummer_NGC362.png

(C1) and (C2) create animations of the orbital evolution of globular clusters
They search for the data file of the orbit as well as the end position of the streams
The end position of the streams stores the escape time
The center of the globular cluster is shown with a red dot in the animations
The size of the red dot changes based on the number of stars that escape

(C3) Creates a nice histogram of the escaped stars, with time, radial, azimuthal and cylindrical radial coordinates

(C4) executes (C1) taking the system argument of the cluster name 