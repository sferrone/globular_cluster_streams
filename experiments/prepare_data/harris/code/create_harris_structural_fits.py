import pandas as pd
from astropy.io import fits
import numpy as np

# inputs
my_frame = pd.read_csv("../inputs/structural_harris_2010.txt", delimiter=';')
my_headers = pd.read_csv("../inputs/structural_harris_2010_header.txt", delimiter=';')
positional_data = fits.open("../outputs/harris_positional_data.fits")
outname = '../outputs/harris_structural_data.fits'
# def variables
outname = '../outputs/harris_structural_data.fits'
my_columns = []
r_c_parsec = np.zeros(positional_data[1].data['ID'].shape)
r_h_l_parsec = np.zeros(positional_data[1].data['ID'].shape)
arcmin_to_radian = np.pi / 60 / 180
counter = 0
# convert arcminutes to parsecs
for x in my_frame.ID:
    #index gymnastics
    positional_index = np.argwhere(positional_data[1].data['ID']==x)[0][0]
    r_c_armin = my_frame.r_c[my_frame['ID']==x].to_numpy()[0]
    r_h_l_arcmin = my_frame.r_h_l[my_frame['ID']==x].to_numpy()[0]
    distance = 1000*positional_data[1].data['R_SUN'][positional_index] # originally in kpc
    r_c_parsec[counter] = arcmin_to_radian*distance*r_c_armin
    r_h_l_parsec[counter] = arcmin_to_radian*distance*r_h_l_arcmin
    counter +=1
# specify data type for each column in the fits tile
my_formats = ["30A", "E", "E", "E", "E", "E", "E" , "E", "E", "E", "E", "E", "E"]
counter=0
# create fits column for each
for x in my_headers['TAG']:
    if x=="r_c":
        my_columns.append(fits.Column(name=x, array=r_c_parsec, format=my_formats[counter], unit="pc"))
    elif x=="r_h_l":
        my_columns.append(fits.Column(name=x, array=r_h_l_parsec, format=my_formats[counter], unit="pc"))
    else:
        my_columns.append(fits.Column(name=x, array=my_frame[x].to_numpy(), format=my_formats[counter], unit=my_headers['units'][counter]))
    counter+=1
# write out table
coldefs = fits.ColDefs(my_columns)
hdu = fits.BinTableHDU.from_columns(coldefs)
hdu.writeto(outname)
positional_data.close()
