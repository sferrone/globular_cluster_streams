import json
import os 

def create_json(outname,time=[-5,5],lon=[-180, 180], lat=[-90,90], D=[0, 200], PMB=[-30, 30], PML_COSB=[-30, 30]):
    """ 
    Purpose: initialize json files for targetting specific glbular clusters 
    All values are defaulted to cover the entire range
    Arguments:
    outname: string
        file name. this should be a stream from ibata 2020. For example "Eridanus" or "NGC3201"
    The rest are coordinates given in galactic coordinates. 
    lat,lon     -> degrees
    D           -> kiloparsec
    PMB         -> Proper motion in b (miliarcseconds / year)
    PML_COSB    -> Proper motion in b (miliarcseconds / year)
    matches     -> leave blank. This will be filled in by a different program
    """
    out_dict = {}
    out_directory = "../outputs/"
    out_dict["time"]        = sorted(time)      # billions of years
    out_dict['LAT']         = sorted(lat)       # degrees
    out_dict['LONG']        = sorted(lon)       # degrees
    out_dict['D']           = sorted(D)         # kpc
    out_dict['PMB']         = sorted(PMB)       # mas / year
    out_dict['PML_COSB']    = sorted(PML_COSB)  # mas / year
    # rewrite the interval file
    fp = open(out_directory+outname+".json", 'w')
    fp.seek(0)
    json.dump(out_dict, fp)
    fp.close()       

if __name__=="__main__":
    """
    create an example output
    """
    lon = [34, 53]
    lat = [26, 75]
    D   = [3, 10]
    PMB = [-9, 2]
    PML_COSB = [-9, 0]
    time = [-1, 1]
    outname = "svol"
    create_json(outname, lat=lat, lon=lon, D=D, PMB=PMB, PML_COSB=PML_COSB,time=time)