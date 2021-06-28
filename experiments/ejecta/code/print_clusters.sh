chmod +x animate_fugitives_rz.py
cat list_of_clusters.txt | while read line
do 
python3 animate_fugitives_rz.py echo $line
done 