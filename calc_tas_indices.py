"""
Filename:    calc_tas_indices.py
Author:      Mitchell Black, mitchell.black@bom.gov.au
Description: Calculate temperature related indices from daily temperature fields 
"""

# Import general Python modules
import xarray as xr
import pandas as pd
import numpy as np  
from datetime import date
from datetime import datetime
import calendar
import seaborn as sns
import time
import sys, os
import glob
import argparse
import git
import cmdline_provenance as cmdprov
import dask.diagnostics
import xesmf as xe
from itertools import chain

import xclim 

# Import my modules

if not os.path.isdir('./modules'):
    sys.path.append('/g/data/mn51/users/mtb563/toolbox/modules')

import utils

# Dask configurations

dask.config.set({"array.slicing.split_large_chunks": True}) 

# Define functions

def get_tasmax(inargs,Ystart,Yend):
    """Read in daily maximum temperature
    Args:
        inargs (class): class object with input arguments from command line
        Y (int): identify year to calculate
    Returns:
        DataArray: data array containing tasmax per sampling constraints.
    """

    fnames = [inargs.tasmax_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)) for Y in range(Ystart,Yend+1)]
    fnames = list(set(fnames))
    fnames.sort()    
    tasmax = utils.generalio.read_data(infiles=fnames,var=inargs.tasmax_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Ystart}-01-01',f'{Yend}-{inargs.yearend}'],output_units='degC')
    
    utils.timeseries.check_correct_ntimesteps(tasmax,sdate=f'{Ystart}-01-01',edate=f'{Yend}-{inargs.yearend}',freq='D')
    
    return tasmax

def get_tasmin(inargs,Ystart,Yend):
    """Read in daily minimum temperature
    Args:
        inargs (class): class object with input arguments from command line
        Y (int): identify year to calculate
    Returns:
        DataArray: data array containing tasmin per sampling constraints.
    """

    fnames = [inargs.tasmin_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)) for Y in range(Ystart,Yend+1)]
    fnames = list(set(fnames))
    fnames.sort()    
    tasmin = utils.generalio.read_data(infiles=fnames,var=inargs.tasmin_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Ystart}-01-01',f'{Yend}-{inargs.yearend}'],output_units='degC')
    
    utils.timeseries.check_correct_ntimesteps(tasmin,sdate=f'{Ystart}-01-01',edate=f'{Yend}-{inargs.yearend}',freq='D')
    
    return tasmin


def model(inargs):
    """Return string summarising model details"""

    model = '_'.join([inargs.driving_model, inargs.pathway, inargs.downscaling_model, inargs.bias_correction_method ])
    
    if inargs.regrid_target_grid:
        model=model+"_regridded"
    
    return model

def get_fname(index,tperiod):
    """Define name for output files created in this program"""
    assert isinstance(index,str)
    assert isinstance(tperiod,str)
    return args.ofile_drs.replace("INDEX",index).replace("TPERIOD",tperiod)

def global_attrs(inargs):
    """Return dictionary of global attributes to add to output file"""
    
    return {
            'description': "Heat indices calculated for the Australian Climate Service",
            'driving_model': inargs.driving_model,
            'downscaling_model': inargs.downscaling_model,
            'pathway': inargs.pathway,
            'bias_correction_method': inargs.bias_correction_method,
            'contact': "Mitchell Black (mitchell.black@bom.gov.au)",
            'code': "https://github.com/Ausutils.climateateService/hazards-heat"
            }

def main(inargs):
    """Calculate the number of days per year when temperatures meet threshold conditions"""

    dask.diagnostics.ProgressBar().register()

    if inargs.calendar == 'standard':
        inargs.yearend = '12-31'
    elif inargs.calendar == '360_day':
        inargs.yearend = '12-30'

    if inargs.ofile_drs:
        if not all(s in inargs.ofile_drs for s in ['INDEX','TPERIOD']):
            raise ValueError("user defined argument ofile_drs must contain strings 'INDEX' and 'TPERIOD'") 
    else:
        inargs.ofile_drs = f'INDEX_{model(inargs)}_day_TPERIOD.nc'

    if not os.path.exists(get_fname(index=inargs.index,tperiod=f"{inargs.StartYr}0101-{inargs.EndYr}{inargs.yearend.replace('-','')}").replace('day','annual')):
        for Y in range(inargs.StartYr,inargs.EndYr+1):
            if not os.path.exists(get_fname(index=inargs.index,tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}').replace('day','annual')):

                if inargs.index == 'TXm':
                    index = xclim.indices.tx_mean(get_tasmax(inargs,Y,Y),freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':'annual mean daily maximum temperature','cell_methods':'time: mean (interval: 1Y)'})
            
                elif inargs.index == 'TXx':
                    index = xclim.indices.tx_max(get_tasmax(inargs,Y,Y),freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':'annual maximum daily maximum temperature','cell_methods':'time: maximum (interval: 1Y)'})
 
                elif inargs.index in ['TXge35','TXge40','TXge45','TXge50']:
                    index = xclim.indices.tx_days_above(get_tasmax(inargs,Y,Y), thresh=f'{float(inargs.index[4:6])} degC', freq='YS', op='>=')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days greater than or equal to {float(inargs.index[4:6])}degC','cell_methods':'time: count (interval: 1Y)'})
                
                elif inargs.index == 'TX90P':
                    if not os.path.exists(get_fname(index='TX90perc_doy',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}')):
                        tasmax_bp = get_tasmax(inargs,inargs.BPStartYr,inargs.BPEndYr)
                        print(tasmax_bp)
                        tasmax_bp_per = xclim.core.calendar.percentile_doy(tasmax_bp,per=90,window=5).sel(percentiles=90)
                        utils.generalio.save_data(tasmax_bp_per,get_fname(index='TX90perc_doy',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))
                        del(tasmax_bp,tasmax_bp_per)
                    
                    tasmax_per = xr.open_dataset(get_fname(index='TX90perc_doy',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))
                    print(tasmax_per)
                    index = xclim.indices.tx90p(tget_tasmax(inargs,Y,Y),tasmax_per)
                    index = update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days above the doy 90th percentile','cell_methods':'time: count (interval: 1Y)'})
                
                elif inargs.index == 'TNle02':
                    index = xclim.indices.tn_days_below(get_tasmin(inargs,Y,Y), thresh=f'{float(inargs.index[4:6])} degC', freq='YS', op='<=')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days less than or equal to {float(inargs.index[4:6])}degC','cell_methods':'time: count (interval: 1Y)'})
                
                utils.generalio.save_data(index,get_fname(index=inargs.index,tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}').replace('day','annual'),append_global_attrs=global_attrs(inargs))
    
        # Concatenate annual files
        fnames = [get_fname(index=inargs.index,tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}').replace('day','annual') for Y in range(inargs.StartYr,inargs.EndYr+1)]
        fnames.sort()    
        index_cat = utils.generalio.read_data(infiles=fnames,var=inargs.index)
        utils.timeseries.check_correct_ntimesteps(index_cat,sdate=f'{inargs.StartYr}-01-01',edate=f'{inargs.EndYr}-{inargs.yearend}',freq='YS')
        utils.generalio.save_data(index_cat,get_fname(index=inargs.index,tperiod=f"{inargs.StartYr}0101-{inargs.EndYr}{inargs.yearend.replace('-','')}").replace('day','annual'),append_global_attrs=global_attrs(inargs))

        if inargs.tidy_wkdir:
            for f in fnames:
                os.remove(f)
    
    print('Analysis complete!')


if __name__ == '__main__':
    extra_info =""" 
author:
    Mitchell Black, mitchell.black@bom.gov.au
"""
    description = """
    Calculate specified indices from daily maximum temperature fields.    
    """
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--index", type=str, choices=['TXm','TXx','TXge35','TXge40','TXge45','TXge50','TX90P','TNle02'], help="specify the index to be computed")
    parser.add_argument("--tasmax_fpath", type=str, default=None, help="generic path to tasmax files (specify year as {Y} and emission pathway as {pathway}")
    parser.add_argument("--tasmax_varname", type=str, default ='tasmax', help="variable name for tasmax in tasmax_fpath")
    parser.add_argument("--tasmin_fpath", type=str, default=None, help="generic path to tasmin files (specify year as {Y} and emission pathway as {pathway}")
    parser.add_argument("--tasmin_varname", type=str, default ='tasmin', help="variable name for tasmin in tasmin_fpath")
    parser.add_argument("--driving_model", type=str, help="Name of the driving model")
    parser.add_argument("--downscaling_model", type=str, help="Name of the downscaling model")
    parser.add_argument("--bias_correction_method", type=str, choices=['raw','qme','ecdfm','mbcn','mrnbc','qdc','ACS-QME'], help="Name of the bias correction method")
    parser.add_argument("--pathway", type=str, choices=['ssp126','ssp370','rcp45','rcp85','historical'], help="Emission pathway")
    parser.add_argument("--BPStartYr", type=int, default=1985, help="Start of EHF base period YYYY")
    parser.add_argument("--BPEndYr", type=int, default=2014, help="End of EHF base period YYYY")
    parser.add_argument("--lon_bounds", type=float, default=None, nargs='*', help="Longitude: single value for nearest point or two values for bounds")
    parser.add_argument("--lat_bounds", type=float, default=None, nargs='*', help="Latitude: single value for nearest point or two values for bounds")
    parser.add_argument("--StartYr", type=int, help="Calculate EHF from this year")
    parser.add_argument("--EndYr", type=int, help="Calculate EHF to this year")
    parser.add_argument("--calendar", type=str, choices=['standard','360_day'], default='standard', help="Calendar type for input data")
    parser.add_argument("--ofile_drs",type=str,default=False,help="Define drs for output files. Must contain INDEX and TPERIOD (replaced by script). Default: INDEX_<driving-model>_<pathway>_<downscaling-model>_<bias-correction_method>_TPERIOD.nc")
    parser.add_argument("--tidy_wkdir",type=bool,default=False,help="Remove intermediate working files")

    args = parser.parse_args()
    main(args)    
