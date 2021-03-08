import pandas as pd
from astropy.io import fits
import numpy as np
# set path data
infile = '../inputs/structural_baumgardt_data.txt'
outfile = 'baumgardt_structural_data.fits'
outpath = "../outputs/"
# import data
my_baumgardt = pd.read_csv(infile, delimiter=';', skiprows=[1])
my_units = pd.read_csv(infile, delimiter=';', skiprows=0, nrows=1)
# convert to fits format
my_columns = []
for x in my_baumgardt:
    if x == "ID":
        my_columns.append(fits.Column(name=x, array=my_baumgardt[x].to_numpy(), format='30A', unit=my_units[x][0]))
    else:
        my_columns.append(fits.Column(name=x, array=my_baumgardt[x].to_numpy(), format='E', unit=my_units[x][0]))
# write the file
coldefs = fits.ColDefs(my_columns)
hdu = fits.BinTableHDU.from_columns(coldefs)
hdu.writeto(outpath+outfile)