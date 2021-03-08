SCOPE:
Create a standard file format for out globular data

Structural data of globular clusters was tajen from (S0) and inserted into (I0). This text file, (I0), was then cleaned by hand in visual studio code. Nan was added to cells with missing data, and regular expressions were used to deliminate each element with a semicolon. The code (C0) then ingested (I0) and put this data into a fits format in (O0).

(C1) is scrap


Source
(S0) https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html
INPUTS/INTERMEDIATE
(I0) inputs/
CODE 
(C0) code/create_baumgardt_fits.py
(C1) code/scrap.ipynb
OUTPUTS
(O0) outputs/baumgardt_structural_data.fits