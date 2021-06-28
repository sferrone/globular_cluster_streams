#!/bin/bash
cat list_of_clusters.txt | while read line
do 
    python3 galactic_plot_individual.py $line galactic
    python3 galactic_plot_individual.py $line equatorial
    python3 galactic_plot_individual.py $line galactocentric
done