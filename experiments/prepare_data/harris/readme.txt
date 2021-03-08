SCOPE:
Text structural and positional of globular clusters from the Harris 2012 catalog and store them as fits files. 

Inputs are the text data taken from (S0) and cleaned with visual studio codes to create (I0) and (I1). 
The positional data was stored in (I0) and the header data was stored in (I1). 
The unix script (C0) elimitaed spaces and prepared lists to be interpreted by python and created (I2)
The python file (C1) ingested the position and header data from (I1) and (I2) and created the fits file (O0)

The structural data in files (I3) was also taken from (S0). It was mostly cleaned by hand, space were removed and semi colons were added as delimiters with this unix command [cat structural_harris_2010.txt | sed 's/ */./g']. (I4) was prepared by hand. These two files were then ingested by (C2) and the output (01) was created.
(C3) is entirely scrap work.

Source
(S0) https://www.physics.mcmaster.ca/~harris/mwgc.dat
INPUTS/INTERMEDIATE
(I0) inputs/positional_data_harris_2010.txt 
(I1) inputs/positional_data_harris_2010_header
(I2) inputs/clean_positional_data_harris.txt
(I3) inputs/structural_harris_2010_header.txt
(I4) inputs/structural_harris_2010.txt
CODE 
(C0) code/make_nice_coordinates.sh
(C1) code/create_harris_positional_fits.py
(C2) code/create_harris_structural_fits.py
(C3) code/Untitled.ipynb
OUTPUTS
(O0) outputs/harris_positional_data.fits
(01) outputs/harris_structural_data.fits