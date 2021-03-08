import pandas as pd
from astropy.io import fits
import numpy as np

# open data
my_frame = pd.read_csv("../inputs/clean_positional_data_harris.txt", delimiter=';')
my_headers = pd.read_csv("../inputs/positional_data_harris_2010_header.txt", delimiter=';')
# define cariables 
my_columns = []
my_formats = ["30A", "3E", "3E", "E", "E", "E", "E" , "E", "E", "E"]
counter=0
temp_array = np.zeros((157,3))
n_data_points = my_frame.shape
# Convert the RA and dec to arrays
def interpret_string(my_string):
    list_string = my_string[1:-2].split(',')
    return np.array([float(x) for x in list_string])
my_frame['DEC']=my_frame['DEC'].apply(interpret_string)
my_frame['RA'] = my_frame['RA'].apply(interpret_string)

# Put data from pandas into fits format
for x in my_headers['TAG']:
    if x=="RA":
        for i in range(n_data_points[0]):
            temp_array[i,:] = my_frame[x][i]        
        my_columns.append(fits.Column(name=x, array=temp_array, format=my_formats[counter], dim='(1,3)', unit=my_headers['units'][counter]))
    elif x=="DEC":
        for i in range(n_data_points[0]):
            temp_array[i,:] = my_frame[x][i]        
        my_columns.append(fits.Column(name=x, array=temp_array, format=my_formats[counter], dim='(1,3)', unit=my_headers['units'][counter]))
    else:
        my_columns.append(fits.Column(name=x, array=my_frame[x].to_numpy(), format=my_formats[counter], unit=my_headers['units'][counter]))
    counter+=1
# Create the file and then close the fits header
t = fits.BinTableHDU.from_columns(my_columns)
t.writeto('../outputs/harris_positional_data.fits')