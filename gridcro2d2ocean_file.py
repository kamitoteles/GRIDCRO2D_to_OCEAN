
#+---------------------------------------------------------------------
#+ Python Script that creates CMAQ OCEAN files form MCIP GRICRO2d file
#+ Check the 'CHANGE' comments to define your directories for the run
#+ Author: Camilo Moreno
#+ Email = cama9709@gmail.com
#+---------------------------------------------------------------------
#%%
import numpy as np
import glob 
from netCDF4 import Dataset
from datetime import datetime, date, timedelta
from os import listdir, scandir, getcwd

def create_forced_var(ds_cro):
    """Create array where cells that wanted to be forced are equal to 1
    and the rest are equal to 0 for multiple locations.

    Keyword arguments:
    ds_cro -- Dataset of the latlon netCDF file of the grid
    dic_coords -- dictionary of latlon coordinates for lower left and 
                  upper right corners of each specified kay-value location
    """
    open_var = ds_cro.variables['LWMASK'][:][:][:][:]
    open_var[open_var > 0] = 2
    open_var[open_var == 0] = 1
    open_var[open_var > 1] = 0

    return open_var

def create_ncfile(save_dir, spc_name, open_var, ds_ocean):
    """Create Final NETCDF file.

    Keyword arguments:
    save_dir -- string of the location for saving the netCDF files
    ds_cro -- Dataset of the latlon netCDF file of the grid
    spc_name -- string of the forced species name
    open_var -- array containing the forced values by location for the species
    ds_ocean -- conc file dataset
    """
    hr = 0
    num_vars = 1
    lays = 1
    ltime = 1
    cols = len(ds_ocean.dimensions['COL'])
    rows = len(ds_ocean.dimensions['ROW'])
    datetimes = len(ds_ocean.dimensions['DATE-TIME'])

    #* Create new netCDF
    new_cmaq_file = f'{save_dir}/ocean_file_2018_NSthAm_CROS.ncf'
    ds_new_cmaq = Dataset(new_cmaq_file, open = True, mode = 'w', format=  "NETCDF3_64BIT")

    #* Create dimenssions
    TSTEP = ds_new_cmaq.createDimension("TSTEP", None)
    DATE_TIME = ds_new_cmaq.createDimension("DATE-TIME", datetimes)
    LAY = ds_new_cmaq.createDimension("LAY", lays)
    VAR = ds_new_cmaq.createDimension("VAR", num_vars)
    ROW = ds_new_cmaq.createDimension("ROW", rows)
    COL = ds_new_cmaq.createDimension("COL", cols)

    ds_new_cmaq.sync()

    #* Creatae attributes
    attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", 
           "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", 
           "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", 
           "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
    
    for attr in attrs:
        if hasattr(ds_ocean, attr): 
            attrVal = getattr(ds_ocean, attr)
            setattr(ds_new_cmaq, attr, attrVal)

    varlist = spc_name + ' '*(16 - len(spc_name))
    cmaq_attrs = {'NLAYS': np.int32(lays),
                  'NVARS': np.int32(num_vars),
                  'UPNAM': "OCEAN_FILE",
                  'VAR-LIST': varlist,
                  'FILEDESC': "Created by Camilo Moreno"
                  }

    for attr in cmaq_attrs:
        ds_new_cmaq.setncattr(attr, cmaq_attrs[attr])

    #* Create variables
    tflag = ds_new_cmaq.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    fill_attrs(ds_ocean, tflag, 'TFLAG')

    var_temp = ds_new_cmaq.createVariable(spc_name,"f4",("TSTEP", "LAY", "ROW", "COL"))
    fill_attrs(ds_ocean, var_temp, 'OPEN')

    #* Fill variables
    ds_new_cmaq.variables[spc_name][:, :, :] = open_var
    dattim = np.squeeze(ds_ocean.variables['TFLAG'][:][:])
    tflag[:] = np.zeros([ltime,lays,datetimes])
    tflag[:] = np.reshape(dattim,[ltime,lays,datetimes])  
    
    #* Close new netcdf file
    ds_new_cmaq.close()

    print(f" OCEAN file DONE")

def fill_attrs(ds_ocean, nc_var, var_name):
    """Fill atribute values for variables as appear in conc file.

    Keyword arguments:
    ds_ocean -- conc file dataset
    nc_var -- netCDF variable of the file
    var_name -- variable name
    """
    varattrs=["long_name","units","var_desc"]
    for varattr in varattrs:
        if hasattr(ds_ocean.variables[var_name], varattr): 
            varattrVal = getattr(ds_ocean.variables[var_name], varattr)
            setattr(nc_var, varattr, varattrVal)

def get_concfiels(ocean_file):
    """Get all files that begins whith 'CONC' from the given directory.

    Keyword arguments:
    ocean_file -- string of the directory whre conc day files
    """
    all_files = [f for f in glob.glob(f'{ocean_file}/CONC.*')]
    all_files.sort()
    return all_files

# %%
if __name__ == "__main__":
    #CHANGE: save_dir: directory path where the forced files will be saved
    #        latlon_file: path of the latlon netCDF of the run grid
    #        ocean_file: path where CONC files are located
    #        spc_name: cbo5_aero5_aq name of the one species to be forced
    #        dic_coord: dictionary containing a list of latlon coords of lower left 
    #                   and upper right corners for each key-value location.
    #                   Format fo dic_coord must be:
    #                   dic_coord = {
    #                               'name_of_location_1': [lower_left_lat_1, lower_left_lon_1, upper_left_lat_1, upper_left_lon_1]
    #                               'name_of_location_2': [lower_left_lat_2, lower_left_lon_2, upper_left_lat_2, upper_left_lon_2]
    #                                 }
    save_dir = '/Volumes/Avispa/OCEAN'
    cro_file = '/Users/camilo/Downloads/ocean/GRIDCRO2D_20180208.nc'
    ocean_file = '/Users/camilo/Downloads/ocean_file_2018_NSthAm_CROS.ncf'
    spc_name = 'OPEN'
    
    ds_cro = Dataset(cro_file, mode = 'r',  open = True)
    ds_ocean = Dataset(ocean_file, mode = 'r',  open = True)

    open_var = create_forced_var(ds_cro)
    create_ncfile(save_dir, spc_name, open_var, ds_ocean)

    ds_ocean.close()
    ds_cro.close()
    
# %%
