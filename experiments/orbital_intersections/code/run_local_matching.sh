#!/bin/bash
cat subset_of_streams.txt | while read line
do 
    python3 intersection_json.py $line 
done