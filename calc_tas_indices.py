"""
Filename:    calc_tas_indices.py
Authors:     Mitchell Black, mitchell.black@bom.gov.au
             Cassandra Rogers, cassandra.rogers@bom.gov.au
Description: Calculate dry-/wet-bulb temperature related indices from daily dry-/wet-bulb temperature fields 
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
        Ystart (int): identify first year to read
        Yend (int): identify last year to read
    Returns:
        DataArray: data array containing tasmax per sampling constraints.
    """

    fnames = [inargs.tasmax_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)) for Y in range(Ystart,Yend+1)]
    fnames = list(set(fnames))
    fnames.sort()    
    tasmax = utils.generalio.read_data(infiles=fnames,var=inargs.tasmax_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Ystart}-01-01',f'{Yend}-{inargs.yearend}'],output_units='degC')

    utils.timeseries.check_correct_ntimesteps(tasmax,sdate=f'{Ystart}-01-01',edate=f'{Yend}-{inargs.yearend}',freq='D')
    
    return tasmax

def get_twisomax(inargs,Ystart,Yend):
    """Read in daily maximum wet-bulb temperature
    Args:
        inargs (class): class object with input arguments from command line
        Ystart (int): identify first year to read
        Yend (int): identify last year to read
    Returns:
        DataArray: data array containing twisomax per sampling constraints.
    """

    fnames = [inargs.twisomax_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)) for Y in range(Ystart,Yend+1)]
    fnames = list(set(fnames))
    fnames.sort()
    twisomax = utils.generalio.read_data(infiles=fnames,var=inargs.twisomax_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Ystart}-01-01',f'{Yend}-{inargs.yearend}'],output_units='degC')

    utils.timeseries.check_correct_ntimesteps(twisomax,sdate=f'{Ystart}-01-01',edate=f'{Yend}-{inargs.yearend}',freq='D')

    return twisomax

def get_tasmin(inargs,Ystart,Yend):
    """Read in daily minimum temperature
    Args:
        inargs (class): class object with input arguments from command line
        Ystart (int): identify first year to read
        Yend (int): identify last year to read
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
    """Calculate the specified index"""

    dask.diagnostics.ProgressBar().register()

    sample_file = inargs.tasmax_fpath.format(Y=inargs.StartYr,pathway=utils.climate.emission_pathway(inargs.StartYr,inargs.pathway))
    nc_calendar = utils.timeseries.get_netcdf_calendar(sample_file)
    
    if nc_calendar == '360_day':
        inargs.yearend = '12-30'
    else:
        inargs.yearend = '12-31'
 
    if inargs.ofile_drs:
        if not all(s in inargs.ofile_drs for s in ['INDEX','TPERIOD']):
            raise ValueError("user defined argument ofile_drs must contain strings 'INDEX' and 'TPERIOD'") 
    else:
        inargs.ofile_drs = f'INDEX_{model(inargs)}_day_TPERIOD.nc'

    if not os.path.exists(get_fname(index=inargs.index,tperiod=f"{inargs.StartYr}0101-{inargs.EndYr}{inargs.yearend.replace('-','')}").replace('day','annual')):
        for Y in range(inargs.StartYr,inargs.EndYr+1):
            if not os.path.exists(get_fname(index=inargs.index,tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}').replace('day','annual')):

                if inargs.index in ['TXm','TXx','TGm','TXge35','TXge40','TXge45','TXge50','TX90P','TX_90P']:
                    max_data = get_tasmax(inargs,Y,Y)
                    max_name = 'temperature'
                elif inargs.index in ['TwXm','TwXx','TwXge25','TwXge27','TwXge29','TwXge31','TwX90P','TwX_90P']:
                    max_data = get_twisomax(inargs,Y,Y)
                    max_name = 'wet-bulb temperature'

                if inargs.index in ['TXm','TwXm']:
                    index = xclim.indices.tx_mean(max_data,freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':f'annual mean daily maximum {max_name}','cell_methods':'time: mean (interval: 1Y)'})
            
                elif inargs.index in ['TXx','TwXx']:
                    index = xclim.indices.tx_max(max_data,freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':f'annual maximum daily maximum {max_name}','cell_methods':'time: maximum (interval: 1Y)'})
                
                elif inargs.index == 'TNm':
                    index = xclim.indices.tn_mean(get_tasmin(inargs,Y,Y),freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':'annual mean daily minimum temperature','cell_methods':'time: mean (interval: 1Y)'})
            
                elif inargs.index == 'TNn':
                    index = xclim.indices.tn_min(get_tasmin(inargs,Y,Y),freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':'annual minimum daily minimum temperature','cell_methods':'time: minimum (interval: 1Y)'})
 
                elif inargs.index == 'TGm':
                    tas = xclim.indices.tas(get_tasmin(inargs,Y,Y),max_data)
                    index = xclim.indices.tg_mean(tas,freq='YS')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':'annual mean daily average temperature','cell_methods':'time: mean (interval: 1Y)'})
                
                elif inargs.index in ['TXge35','TXge40','TXge45','TXge50','TwXge25','TwXge27','TwXge29','TwXge31']:
                    deg_text = inargs.index.replace('TwXge','').replace('TXge','')
                    index = xclim.indices.tx_days_above(max_data, thresh=f'{float(deg_text)} degC', freq='YS', op='>=')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days greater than or equal to {float(deg_text)}degC','cell_methods':'time: count (interval: 1Y)'})
                
                elif inargs.index in ['TX90P','TwX90P']:
                    if inargs.index == 'TX90P':
                        bp_index_name = 'TX90perc_doy'
                        index_varname = inargs.tasmax_varname
                    elif inargs.index == 'TwX90P':
                        bp_index_name = 'TwX90perc_doy'
                        index_varname = inargs.twisomax_varname
                    if not os.path.exists(get_fname(index=bp_index_name,tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}')):
                        print('Creating base period file')
                        if inargs.index == 'TX90P':
                            max_bp = get_tasmax(inargs,inargs.BPStartYr,inargs.BPEndYr)
                        elif inargs.index == 'TwX90P':
                            max_bp = get_twisomax(inargs,inargs.BPStartYr,inargs.BPEndYr)
                        print(max_bp)
                        max_bp_per = xclim.core.calendar.percentile_doy(max_bp,per=90,window=5).sel(percentiles=90)
                        utils.generalio.save_data(max_bp_per,get_fname(index=bp_index_name,tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))
                        del(max_bp,max_bp_per)
                        
                    max_per = xr.open_dataset(get_fname(index=bp_index_name,tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))
                    print(max_per)
                    index = xclim.indices.tx90p(max_data,max_per.per)
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days above the doy 90th percentile (base period={inargs.BPStartYr}-{inargs.BPEndYr}, window=5)','cell_methods':'time: count (interval: 1Y)'})
                
                elif inargs.index == 'TNle02':
                    index = xclim.indices.tn_days_below(get_tasmin(inargs,Y,Y), thresh=f'{float(inargs.index[4:6])} degC', freq='YS', op='<=')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'1',\
                            'long_name':f'days less than or equal to {float(inargs.index[4:6])}degC','cell_methods':'time: count (interval: 1Y)'})

                elif inargs.index in ['TX_90P','TwX_90P']:
                    set_percentile = inargs.index.replace('TX_','').replace('TwX_','').replace('P','')
                    index = max_data.resample(time='YS').quantile(float(set_percentile)/100,dim='time',skipna=True,keep_attrs=True,method='midpoint')
                    index = utils.generalio.update_attrs(index,{'name':inargs.index,'units':'degC',\
                            'long_name':f'{float(set_percentile)}th percentile of {max_name}','cell_methods':f'time: {float(set_percentile)}th percentile (interval: 1Y)'})
                
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
authors:
    Mitchell Black, mitchell.black@bom.gov.au
    Cassandra Rogers, cassandra.rogers@bom.gov.au
"""
    description = """
    Calculate specified indices from daily maximum/minimum temperature or daily maximum wet-bulb temperature fields.    
    """
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--index", type=str, choices=['TNm','TNn','TXm','TXx','TGm','TXge35','TXge40','TXge45','TXge50','TX90P','TX_90P','TNle02','TwXm','TwXx','TwXge25','TwXge27','TwXge29','TwXge31','TwX90P','TwX_90P'], help="specify the index to be computed")
    parser.add_argument("--tasmax_fpath", type=str, default=None, help="generic path to tasmax files (specify year as {Y} and emission pathway as {pathway}")
    parser.add_argument("--tasmax_varname", type=str, default='tasmax', help="variable name for tasmax in tasmax_fpath")
    parser.add_argument("--twisomax_fpath", type=str, default=None, help="generic path to twisomax files (specify year as {Y} and emission pathway as {pathway}")
    parser.add_argument("--twisomax_varname", type=str, default='twisomax', help="variable name for twisomax in twisomax_fpath")
    parser.add_argument("--tasmin_fpath", type=str, default=None, help="generic path to tasmin files (specify year as {Y} and emission pathway as {pathway}")
    parser.add_argument("--tasmin_varname", type=str, default='tasmin', help="variable name for tasmin in tasmin_fpath")
    parser.add_argument("--driving_model", type=str, help="Name of the driving model")
    parser.add_argument("--downscaling_model", type=str, help="Name of the downscaling model")
    parser.add_argument("--bias_correction_method", type=str, choices=['raw','qme','ecdfm','mbcn','mrnbc','qdc','ACS-QME','ACS-MRNBC','QDC'], help="Name of the bias correction method")
    parser.add_argument("--pathway", type=str, choices=['ssp126','ssp370','rcp45','rcp85','historical'], help="Emission pathway")
    parser.add_argument("--BPStartYr", type=int, default=1985, help="Start of index base period YYYY")
    parser.add_argument("--BPEndYr", type=int, default=2014, help="End of index base period YYYY")
    parser.add_argument("--lon_bounds", type=float, default=None, nargs='*', help="Longitude: single value for nearest point or two values for bounds")
    parser.add_argument("--lat_bounds", type=float, default=None, nargs='*', help="Latitude: single value for nearest point or two values for bounds")
    parser.add_argument("--StartYr", type=int, help="Calculate index from this year")
    parser.add_argument("--EndYr", type=int, help="Calculate index to this year")
    parser.add_argument("--ofile_drs",type=str,default=False,help="Define drs for output files. Must contain INDEX and TPERIOD (replaced by script). Default: INDEX_<driving-model>_<pathway>_<downscaling-model>_<bias-correction_method>_TPERIOD.nc")
    parser.add_argument("--tidy_wkdir",type=bool,default=False,help="Remove intermediate working files")

    args = parser.parse_args()
    main(args)    
