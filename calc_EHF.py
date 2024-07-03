"""
Filename:    calc_EHF.py
Author:      Mitchell Black, mitchell.black@bom.gov.au
Description: Calculate the Excess Heat Factor 
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
# Import my modules

if not os.path.isdir('./modules'):
    sys.path.append('/g/data/mn51/users/mtb563/toolbox/modules')

import utils

# Dask configurations

dask.config.set({"array.slicing.split_large_chunks": True}) 

# Define functions

def calc_EHF(T,T95):
    """
    Calculate the Excess Heat Factor
    Args:
        T (xarray dataarray): timeseries of daily mean temperature
        T95 (xarray dataarray): 95th percentile of daily mean temperature
    Returns:
        xarray dataarray: daily Excess Heat Factor   
    """

    EHI_accl = xr.DataArray(data=np.nan,coords=T.coords,dims=T.dims)
    EHI_sig = xr.DataArray(data=np.nan,coords=T.coords,dims=T.dims)
    
    for i in range(32,len(T),1):
        EHI_accl[i,...] = (T[i-2:i+1,...].sum(dim='time') / 3) - (T[i-32:i-2,...].sum(dim='time') / 30)  
        EHI_sig[i,...] = (T[i-2:i+1,...].sum(dim='time') / 3) - T95
    
    EHF = EHI_sig * EHI_accl.where(EHI_accl>1,1)

    EHI_accl = utils.generalio.update_attrs(EHI_accl,{'name':'EHI_accl','long_name':'Excess Heat Acclimatisation Index','standard_name':'excess_heat_acclimatisation_index','units':'K'})
    EHI_sig = utils.generalio.update_attrs(EHI_sig,{'name':'EHI_sig','long_name':'Excess Heat Significance Index','standard_name':'excess_heat_signficance_index','units':'K'})
   
    EHF = utils.generalio.update_attrs(EHF,{'name':'EHF','long_name':'Excess Heat Factor','standard_name':'excess_heat_factor','units':'K**2'})
    ds = xr.merge([EHF,EHI_accl,EHI_sig])
    ds.attrs = {}
    return ds

def calc_T95(inargs):
    """Calculate the 95th percentile of daily mean temperature over specified base period.
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
    
    Returns:
        None : saves xarray dataset to NetCDF file
    
    """
    
    print(f'Calculating 95th percentile of daily mean temperature, {inargs.BPStartYr} - {inargs.BPEndYr}')
    
    fnames = [get_fname(var='tasmean',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}') for Y in range(inargs.BPStartYr,inargs.BPEndYr+1)]
    fnames.sort()
    ds_mf = xr.open_mfdataset(fnames, chunks={'time': '100MB'}, concat_dim="time", combine="nested",data_vars='minimal', coords='minimal', compat='override', parallel=True)
    Tb = ds_mf['tasmean'].chunk(-1).sortby('time')
    
    utils.timeseries.check_correct_ntimesteps(Tb,f'{inargs.BPStartYr}-01-01',f'{inargs.BPEndYr}-{inargs.yearend}',freq='D')

    T95 = Tb.quantile(q=0.95,dim='time',keep_attrs=True)
    utils.generalio.save_data(T95,get_fname(var='T95',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'),append_global_attrs=global_attrs(inargs))
    
    del(Tb,T95)
    
def calc_EHF85(inargs):
    """Calculate the 85th percentile of daily EHF values over specified base period.
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
    
    Returns:
        None : saves xarray dataset to NetCDF file
    
    """

    print(f'Calculating EHF 85th percentile, {inargs.BPStartYr} - {inargs.BPEndYr}')

    fnames = [get_fname(var='EHF',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}') for Y in range(inargs.BPStartYr,inargs.BPEndYr+1)]
    fnames.sort()
    ds_mf = xr.open_mfdataset(fnames, chunks={'time': '100MB'}, concat_dim="time", combine="nested",data_vars='minimal', coords='minimal', compat='override', parallel=True)  
    EHF = ds_mf['EHF'].chunk(-1)
    
    utils.timeseries.check_correct_ntimesteps(EHF,f'{inargs.BPStartYr}-01-01',f'{inargs.BPEndYr}-{inargs.yearend}',freq='D')

    EHF85 = EHF.where(EHF>0).quantile(q=0.85,dim='time')
    utils.generalio.save_data(EHF85,get_fname(var='EHF85',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'),append_global_attrs=global_attrs(inargs))
    del(ds_mf,EHF,EHF85)

def calc_EHF_Year(inargs,Y):
    """Calculate the Excess Heat Factor for a specified year.
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
        Y (int): identify year to calculate
    Returns:
        xarray dataset: daily Excess Heat Factor for specified year
    """

    print('Calculating EHF:', Y)
    
    if Y == inargs.StartYr and inargs.ShortStartYr:
        ds = xr.open_dataset(get_fname('tasmean',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}'))
        tasmean = ds['tasmean'].sortby('time').sel(time=slice(f'{Y}-01-01',f'{Y}-{inargs.yearend}')).load()
    else:
        ds = xr.open_mfdataset([get_fname(var='tasmean',tperiod=f'{YY}0101-{YY}{inargs.yearend.replace("-","")}') for YY in [Y-1,Y]],concat_dim='time',combine="nested")
        tasmean = ds['tasmean'].sortby('time').sel(time=slice(f'{Y-1}-11-25',f'{Y}-{inargs.yearend}')).load()
    
    T95 = xr.open_dataset(get_fname(var='T95',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))['tasmean'] 
    EHF = calc_EHF(tasmean,T95)
    EHF = EHF.sel(time=slice(f'{Y}-01-01',f'{Y}-{inargs.yearend}'))
    utils.generalio.save_data(EHF,get_fname(var='EHF',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}'),append_global_attrs=global_attrs(inargs))
    del(ds,tasmean,T95,EHF)

def identify_heatwave_days(EHF, EHF_calc_window = 'inc_prior_days'):
    """Identify all days within the three day periods that contribute to positive EHF values (i.e., heatwaves)
    Args:
        EHF (xarray dataarray): timeseries of EHF values
        EHF_calc_window: speficies if EHF calculated using three day periods including two 'prior' or 'future' days.
    Returns:
        xarray dataarray: timeseries identifying heatwave days (hw = 1., no = 0.)
    """
    event_trigger = xr.where(EHF > 0,1.,0.)
    
    if EHF_calc_window == 'inc_prior_days':
        HW_day = event_trigger.reindex(time=event_trigger.time[::-1]).rolling(time=3,min_periods=1).max()
        HW_day = HW_day.reindex(time=event_trigger.time)
    
    elif EHF_calc_window == 'inc_future_days':
        HW_day = event_trigger.rolling(time=3,min_periods=1).max()
    
    return HW_day

def get_days_EHF_rating(EHFsev,rating):
    """Identify days where the EHF 'rating' meetings the required rating
    Args:
        EHFsev (xarray dataarray): timeseries of gridded EHFsev values (EHFsev=EHF/EHF85)
        rating (str): heatwave classification [Low-Moderate, Severe or Extreme]
    Returns:
        xarray dataarray: timeseries identifying days with EHF meeting desired rating (true = 1.)
    """
    if rating == 'Low-Moderate':
        return xr.where( (EHFsev > 0) & (EHFsev < 1), EHFsev, np.nan)
    if rating == 'Severe':
        return xr.where( (EHFsev >= 1) & (EHFsev < 3), EHFsev, np.nan)
    if rating == 'Extreme':
        return xr.where( EHFsev >= 3, EHFsev, np.nan) 
    
def calc_annual_hw_metrics(inargs,Y):
    """Calculate annual heatwave metrics from gridded EHF data
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
        Y (int): identify year to calculate
    Returns:
        None : saves xarray dataset to NetCDF file
    """
    print('Calculating annual heatwave metrics:', Y)

    if inargs.annual_groupby == 'Jan2Dec':
        Y = [Y] 
        sdate,edate = f'{Y[0]}-01-01', f'{Y[0]}-{inargs.yearend}'
    elif inargs.annual_groupby == 'Jul2Jun':
        Y = [Y,Y+1]
        sdate,edate = f'{Y[0]}-07-01', f'{Y[1]}-06-30'
    
    EHF = utils.generalio.read_data(infiles=[get_fname(var='EHF',tperiod=f'{YY}0101-{YY}{inargs.yearend.replace("-","")}') for YY in Y],var='EHF',\
            lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[sdate,edate])
    
    utils.timeseries.check_correct_ntimesteps(EHF,sdate,edate,freq='D')

    tasmax = utils.generalio.read_data(infiles=list(dict.fromkeys([inargs.tasmax_fpath.format(Y=YY,pathway=utils.climate.emission_pathway(YY,inargs.pathway)) for YY in Y])),\
            var=inargs.tasmax_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[sdate,edate],output_units='degC')
    
    utils.timeseries.check_correct_ntimesteps(tasmax,sdate,edate,freq='D')

    if inargs.regrid_target_grid:
        target_grid = utils.generalio.read_data(infiles=inargs.regrid_target_grid,var=None,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,return_as_DataArray=False)
        tasmax = utils.generalio.regrid(tasmax,target_grid,method=inargs.regrid_method)

    EHF85 = xr.open_dataset(get_fname(var='EHF85',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}'))['EHF']

    resamp_alias = {'Jan2Dec':'YS-JAN', 'Jul2Jun':'YS-JUL'}
    
    HW_day = identify_heatwave_days(EHF,EHF_calc_window = 'inc_prior_days')
    HW_day = HW_day.compute()

    #HWF - frequency of heatwave days in a given year.
    HWF = utils.timeseries.count(HW_day.where(HW_day > 0).to_dataset(name='HWF'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWF = utils.generalio.update_attrs(HWF.HWF,{'units':'1','long_name':'number of days experiencing heatwave conditions'})
    HWF = HWF.compute()

    #HWD - duration of the longest heatwave in a given year.
    HW_dur = utils.timeseries.cumsum_reset_at_value(HW_day,reset_value=0)
    HWD = utils.timeseries.maximum(HW_dur.to_dataset(name='HWD'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWD = utils.generalio.update_attrs(HWD.HWD,{'units':'1','long_name':'duration of longest heatwave (days)'})
    HWD = HWD.compute()    

    #HWN - total number of separate heatwaves in a given year,
    HWN = utils.timeseries.count(HW_dur.where(HW_dur == 1.).to_dataset(name='HWN'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWN = utils.generalio.update_attrs(HWN.HWN,{'units':'1','long_name':'number of separate heatwave events'})
    HWN = HWN.compute()

    #HWAtx - heatwave amplitude - hottest temperature in a given year
    HWAtx = utils.timeseries.maximum(tasmax.where(HW_day == 1.).to_dataset(name='HWAtx'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWAtx = utils.generalio.update_attrs(HWAtx.HWAtx,{'units':'degC','long_name':'maximum temperature of hottest heatwave'})
    HWAtx = HWAtx.compute()

    #HWAsev - heatwave amplitude - highest severity rating in a given year
    EHF_sev = EHF / EHF85
    EHF_sev = EHF_sev.drop_vars('quantile')
    HWAsev = utils.timeseries.maximum(EHF_sev.to_dataset(name='HWAsev'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWAsev = utils.generalio.update_attrs(HWAsev.HWAsev,{'units':'1',\
            'long_name':'Maximum heatwave severity rating: "low-intensity" (0< HWAsev < 1), "severe" (1<= HWAsev < 3), "extreme" (3<= HWAsev)'})
    HWAsev = HWAsev.compute()

    #HW_low - number of days with EHF rating of 'Low-Moderate'
    HW_low = utils.timeseries.count(get_days_EHF_rating(EHF_sev,'Low-Moderate').to_dataset(name='HW_low'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HW_low = utils.generalio.update_attrs(HW_low.HW_low,{'units':'1',\
            'long_name':'Number of three day periods (TDPs) with "Low-Moderate" heatwave conditions'})
    HW_low = HW_low.compute()

    #HW_sev - number of days with EHF rating of 'Severe'
    HW_sev = utils.timeseries.count(get_days_EHF_rating(EHF_sev,'Severe').to_dataset(name='HW_sev'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HW_sev = utils.generalio.update_attrs(HW_sev.HW_sev,{'units':'1',\
            'long_name':'Number of three day periods (TDPs) with "Severe" heatwave conditions'})
    HW_sev = HW_sev.compute()

    #HW_ext - number of days with EHF rating of 'Extreme'
    HW_ext = utils.timeseries.count(get_days_EHF_rating(EHF_sev,'Extreme').to_dataset(name='HW_ext'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HW_ext = utils.generalio.update_attrs(HW_ext.HW_ext,{'units':'1',\
            'long_name':'No. of three day periods (TDPs) with "Extreme" heatwave conditions'})
    HW_ext = HW_ext.compute()

    days_since = xr.DataArray(data=np.ones(HW_day.shape),dims=HW_day.dims,coords=HW_day.coords) * xr.DataArray(data=list(range(len(HW_day.time))),dims=('time'),coords={'time':HW_day.time})
    days_since = days_since.where(HW_day == 1.0)
    days_since = days_since.compute()

    #HWTf - timing of first heatwave in season  
    HWTf = utils.timeseries.minimum(days_since.to_dataset(name='HWTf'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWTf = utils.generalio.update_attrs(HWTf.HWTf,{'units': f'days from 1 {inargs.annual_groupby[0:3]}',\
            'long_name':'heatwave timing - start of the first heatwave'})
    HWTf = HWTf.compute()

    #HWTl - timing of last heatwave in season
    HWTl = utils.timeseries.maximum(days_since.to_dataset(name='HWTl'),resamp={'time':resamp_alias[inargs.annual_groupby]})
    HWTl = utils.generalio.update_attrs(HWTl.HWTl,{'units': f'days from 1 {inargs.annual_groupby[0:3]}',\
            'long_name':'heatwave timing - end of the last heatwave'})
    HWTl = HWTl.compute()

    da_out = [HWF,HWD,HWN,HWAtx,HWAsev,HW_low,HW_sev,HW_ext,HWTf,HWTl]
    for da in da_out:
        da = da.where(~EHF85.isnull())
    ds_out = xr.merge(da_out)
    ds_out.attrs = {}
    ds_out.time.encoding['units'] = "days since 1970-01-01 00:00:00"
    ds_out['time'].attrs["bounds"] = "time_bnds"
    utils.generalio.save_data(ds_out,get_fname(var='EHF-metrics',tperiod=f"{sdate.replace('-','')}-{edate.replace('-','')}_{inargs.annual_groupby}").replace('day','annual'),append_global_attrs=global_attrs(inargs))
    del(HWF,HWD,HWN,HWAtx,HWAsev,HW_low,HW_sev,HW_ext,HWTf,HWTl)

def concat_annual_hw_metrics(inargs):
    """Combine annual heatwave metric files into a single file
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
    Returns:
        None : saves xarray dataset to NetCDF file
    """

    if inargs.annual_groupby == 'Jan2Dec':
        sdate, edate = f'{inargs.StartYr}-01-01', f'{inargs.EndYr}-{inargs.yearend}'
    elif inargs.annual_groupby == 'Jul2Jun':
        sdate, edate = f'{inargs.StartYr}-07-01', f'{inargs.EndYr}-06-30'

    fnames = glob.glob(get_fname(var='EHF-metrics',tperiod=f'*_{inargs.annual_groupby}').replace('day','annual'))
    fnames.sort()
    ds_mf = xr.open_mfdataset(fnames)
    ds_mf = ds_mf.sortby('time').sel(time=slice(sdate,edate))
    utils.generalio.save_data(ds_mf,get_fname(var='EHF-metrics',tperiod=f"{sdate.replace('-','')}-{edate.replace('-', '')}_{inargs.annual_groupby}").replace('day','annual'),append_global_attrs=global_attrs(inargs))
    
    for f in fnames:
        os.remove(f)

def tas_data_exists(inargs,Y):
    """Check if temperature data exists for specified year
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
        Y (int): identify year to calculate
    Returns:
        boolean : True/False depending on data status
    """

    # Consider multiple input files (split per year)
    if '{Y}' in inargs.tasmax_fpath:
        if os.path.exists(inargs.tasmax_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway))):
            return True
        else:
            return False

    # Consider a single input file (with concatenated years)
    if Y in xr.open_dataset(inargs.tasmax_fpath.format(pathway=utils.climate.emission_pathway(Y,inargs.pathway))).time.dt.year:
        return True
    else:
        return False

def calc_tasmean(inargs,Y):
    """Calculate daily mean temperature
    Args:
        inargs (class): class object with input arguments from calc_EHF command line
        Y (int): identify year to calculate
    Returns:
        None : saves xarray dataset to NetCDF file
    """

    if (Y == inargs.StartYr - 1) and (Y not in range(inargs.BPStartYr-1,inargs.BPEndYr+1)):
        if not tas_data_exists(inargs,Y):
            if inargs.annual_groupby == "Jul2Jun":
                inargs.ShortStartYr = True
                return
    
    if not tas_data_exists(inargs,Y):
        raise ValueError(f'Input temperature data does not exist for year {Y}')

    tasmax = utils.generalio.read_data(infiles=inargs.tasmax_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)),\
            var=inargs.tasmax_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Y}-01-01',f'{Y}-{inargs.yearend}'],output_units='degC')
    tasmin = utils.generalio.read_data(infiles=inargs.tasmin_fpath.format(Y=Y,pathway=utils.climate.emission_pathway(Y,inargs.pathway)),\
            var=inargs.tasmin_varname,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,time_bounds=[f'{Y}-01-01',f'{Y}-{inargs.yearend}'],output_units='degC')
    tasmean = (tasmin + tasmax) / 2
    tasmean = utils.generalio.update_attrs(tasmean,{'name':'tasmean','units':'degC','standard_name':'daily_mean_temperature','long_name':'daily_mean_temperature'})
    
    if inargs.regrid_target_grid:
        target_grid = utils.generalio.read_data(infiles=inargs.regrid_target_grid,var=None,lat_bounds=inargs.lat_bounds,lon_bounds=inargs.lon_bounds,return_as_DataArray=False)
        tasmean = utils.generalio.regrid(tasmean,target_grid,method=inargs.regrid_method)

    utils.generalio.save_data(tasmean,get_fname(var='tasmean',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}'),append_global_attrs=global_attrs(inargs))
    del(tasmean,tasmin,tasmax)

def model(inargs):
    """Return string summarising model details"""

    model = '_'.join([inargs.driving_model, inargs.pathway, inargs.downscaling_model, inargs.bias_correction_method ])
    
    if inargs.regrid_target_grid:
        model=model+"_regridded"
    
    return model

def get_fname(var,tperiod):
    """Define name for output files created in this program"""
    assert isinstance(var,str)
    assert isinstance(tperiod,str)
    return args.ofile_drs.replace("VAR",var).replace("TPERIOD",tperiod)

def global_attrs(inargs):
    """Return dictionary of global attributes to add to output file"""
    
    return {
            'description': "Heatwave metrics as defined using the Excess Heat Factor index of Nairn and Fawcett (2015)",
            'driving_model': inargs.driving_model,
            'downscaling_model': inargs.downscaling_model,
            'pathway': inargs.pathway,
            'bias_correction_method': inargs.bias_correction_method,
            'contact': "Mitchell Black (mitchell.black@bom.gov.au)",
            'code': "https://github.com/AusClimateService/hazards-heat"
            }

def main(inargs):
    """Calculate the Excess Heat Factor from gridded model data"""

    dask.diagnostics.ProgressBar().register()

    inargs.ShortStartYr = False
    
    if inargs.calendar == 'standard':
        inargs.yearend = '12-31'
    elif inargs.calendar == '360_day':
        inargs.yearend = '12-30'

    if inargs.ofile_drs:
        if not all(s in inargs.ofile_drs for s in ['VAR','TPERIOD']):
            raise ValueError("user defined argument ofile_drs must contain strings 'VAR' and 'TPERIOD'") 
    else:
        inargs.ofile_drs = f'VAR_{model(inargs)}_day_TPERIOD.nc'
    
    # Compute T95 file
    if not os.path.exists(get_fname(var='T95',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}')):
        for Y in range(inargs.BPStartYr,inargs.BPEndYr+1):
            if not os.path.exists(get_fname(var='tasmean',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}')):
                calc_tasmean(inargs,Y)
        calc_T95(inargs)

    # Compute EHF85 file        
    if not os.path.exists(get_fname(var='EHF85',tperiod=f'{inargs.BPStartYr}0101-{inargs.BPEndYr}{inargs.yearend.replace("-","")}')):
        for Y in range(inargs.BPStartYr,inargs.BPEndYr+1):
            if not os.path.exists(get_fname(var='EHF',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}')):
                for YY in [Y-1,Y]:
                    if not os.path.exists(get_fname(var='tasmean',tperiod=f'{YY}0101-{YY}{inargs.yearend.replace("-","")}')):
                        calc_tasmean(inargs,YY)
                calc_EHF_Year(inargs,Y)
        calc_EHF85(inargs)
    
    # Compute EHF for remaining years
    for Y in range(inargs.StartYr,inargs.EndYr+1):
        if not os.path.exists(get_fname(var='EHF',tperiod=f'{Y}0101-{Y}{inargs.yearend.replace("-","")}')):
            for YY in [Y-1,Y]:
                if not os.path.exists(get_fname(var='tasmean',tperiod=f'{YY}0101-{YY}{inargs.yearend.replace("-","")}')):
                    calc_tasmean(inargs,YY)
            calc_EHF_Year(inargs,Y)
    
    # Compute annual heatwave metrics
    if inargs.annual_groupby == 'Jan2Dec':
        sMMDD, eMMDD = '0101', f'{inargs.yearend.replace("-","")}'
    elif inargs.annual_groupby == 'Jul2Jun':
        sMMDD, eMMDD = '0701', f'0630'

    if not os.path.exists(get_fname(var='EHF-metrics',tperiod=f"{inargs.StartYr}{sMMDD}-{inargs.EndYr}{eMMDD}_{inargs.annual_groupby}").replace('day','annual')):
        for Y in range(inargs.StartYr,inargs.EndYr+(1 if inargs.annual_groupby=="Jan2Dec" else 0)):
            if not os.path.exists(get_fname(var='EHF-metrics',tperiod=f"{Y}{sMMDD}-{Y+1 if inargs.annual_groupby == 'Jul2Jun' else Y}{eMMDD}_{inargs.annual_groupby}").replace('day','annual')):
                calc_annual_hw_metrics(inargs,Y)
        
        concat_annual_hw_metrics(inargs)
        
        if inargs.tidy_wkdir:
            types = ('tasmean*.nc', f'EHF_*.nc')            
            rm_files = []
            for t in types:
                rm_files.extend(glob.glob(t))
            for f in rm_files:
                os.remove(f)
    
    print('Analysis complete!')


if __name__ == '__main__':
    extra_info =""" 
author:
    Mitchell Black, mitchell.black@bom.gov.au
"""
    description = """
    Calculate the Excess Heat Factor from gridded model data. 
    This code is based on Nairn and Fawcett (2015):
        'The Excess Heat Factor: A metric for heatwave intensity and its use in classifying heatwave severity', Int. J. Environ. Res. Public Health 2015, 12 
    Note: this code calculates EHF in hindcast mode (instead of forecast mode used by Nairn and Fawcett).    
        """
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
                                     
    parser.add_argument("--tasmax_fpath", type=str, help="generic path to tasmax files (specify year as {Y} and pathway as {pathway}")
    parser.add_argument("--tasmax_varname", type=str, default ='tasmax', help="variable name for tasmax in tasmax_fpath")
    parser.add_argument("--tasmin_fpath", type=str, help="generic path to tasmin files (specify year as {Y} and pathway as {pathway}")
    parser.add_argument("--tasmin_varname", type=str, default ='tasmin', help="variable name for tasmin in tasmin_fpath")
    parser.add_argument("--driving_model", type=str, help="Name of the driving model")
    parser.add_argument("--downscaling_model", type=str, help="Name of the downscaling model")
    parser.add_argument("--bias_correction_method", type=str, choices=['raw','qme','ecdfm','mbcn','mrnbc','qdc','ACS-QME'], help="Name of the bias correction method")
    parser.add_argument("--pathway", type=str, choices=['ssp126','ssp370','rcp45','rcp85','historical'], help="Emission pathway")
    parser.add_argument("--BPStartYr", type=int, default=1985, help="Start of EHF base period YYYY")
    parser.add_argument("--BPEndYr", type=int, default=2014, help="End of EHF base period YYYY")
    parser.add_argument("--StartYr", type=int, help="Calculate EHF from this year")
    parser.add_argument("--EndYr", type=int, help="Calculate EHF to this year")
    parser.add_argument("--lon_bounds", type=float, default=None, nargs='*', help="Longitude: single value for nearest point or two values for bounds")
    parser.add_argument("--lat_bounds", type=float, default=None, nargs='*', help="Latitude: single value for nearest point or two values for bounds")
    parser.add_argument("--annual_groupby", type=str, choices=['Jan2Dec','Jul2Jun'], help="12 month period (grouping) for calculating annual metrics")
    parser.add_argument("--regrid_target_grid", type=str, default=None, help="Path to dataset containing target horizontal grid")
    parser.add_argument("--regrid_method", type=str, choices=['bilinear','conservative'], help="Regridding method")
    parser.add_argument("--calendar", type=str, choices=['standard','360_day'], default='standard', help="Calendar type for input data")
    parser.add_argument("--ofile_drs",type=str,default=False,help="Define drs for output files. Must contain VAR and TPERIOD (replaced by script). Default: VAR_<driving-model>_<pathway>_<downscaling-model>_<bias-correction_method>_TPERIOD.nc")
    parser.add_argument("--tidy_wkdir",type=bool,default=False,help="Remove intermediate working files")

    args = parser.parse_args()
    main(args)    
