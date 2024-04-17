""" Functions related to operations on time axis """

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

def check_correct_ndays(ds,sdate,edate,inclusive='both'):
    """Check that the dataset/array has the correct number of days.
    Args:
         ds (xarray dataset or dataarray): Dataset with time axis to check
         sdate (str): target start date 'YYYY-MM-DD' 
         edate (str): target end date 'YYYY-MM_DD'
         inclusive (str): default = 'both' (includes sdate and edate in nday count)
     Returns:
         raises AssertionError if len (ds.time) != target_ndays
    """ 
    try:
        calendar = ds.time.to_index().calendar
    except AttributeError:
        calendar = 'standard'
    
    target_ndays = xr.cftime_range(sdate,edate,calendar=calendar,inclusive=inclusive).size
    
    assert ds.time.size == target_ndays, f'dataset does not have sufficient number of days, {ds.time.size} instead of {target_ndays}'

def maximum(ds, resamp, exclude_list=["time_bnds"]):
    """Get maximum from dataset and add time_bounds information.
    Args:
        ds (xarray dataset): Input dataset
        resamp (dict): resampling (e.g., {time:1D}).
    Returns:
        xarray dataset : The resampled xarray dataset
    """
    operation = "maximum"
    ds_out = add_cell_methods(ds, resamp=resamp, operation=operation, exclude_list=exclude_list)
    ds_out = ds_out.resample(resamp).max(keep_attrs=True) if resamp is not None else ds_out.max()
    ds_out = update_time_bnds(ds, ds_out, resamp)
    return ds_out

def minimum(ds, resamp, exclude_list=["time_bnds"]):
    """Get minimum from dataset and add time_bounds information.
    Args:
        ds (xarray dataset): Input dataset
        resamp (dict): resampling (e.g., {time:1D}).
    Returns:
        xarray dataset : The resampled xarray dataset
    """
    operation = "minimum"
    ds_out = add_cell_methods(ds, resamp=resamp, operation=operation, exclude_list=exclude_list)
    ds_out = ds_out.resample(resamp).min(keep_attrs=True) if resamp is not None else ds_out.min()
    ds_out = update_time_bnds(ds, ds_out, resamp)
    return ds_out

def mean(ds, resamp, exclude_list=["time_bnds"]):
    """Get mean from dataset and add time_bounds information.
    Args:
        ds (xarray dataset): Input dataset
        resamp (dict): resampling (e.g., {time:1D}).
    Returns:
        xarray dataset : The resampled xarray dataset
    """
    operation = "mean"
    ds_out = add_cell_methods(ds, resamp=resamp, operation=operation, exclude_list=exclude_list)
    ds_out = ds_out.resample(resamp).mean(keep_attrs=True) if resamp is not None else ds_out.mean()
    ds_out = update_time_bnds(ds, ds_out, resamp)
    return ds_out

def count(ds, resamp, exclude_list=["time_bnds"]):
    """Get count from dataset and add time_bounds information.
    Args:
        ds (xarray dataset): Input dataset
        resamp (dict): resampling (e.g., {time:1D}).
    Returns:
        xarray dataset : The resampled xarray dataset
    """
    operation = "count"
    ds_out = add_cell_methods(ds, resamp=resamp, operation=operation, exclude_list=exclude_list)
    ds_out = ds_out.resample(resamp).count(keep_attrs=True) if resamp is not None else ds_out.count()
    ds_out = update_time_bnds(ds, ds_out, resamp)
    return ds_out

def ndays_since_MMDD(date,ref_month,ref_day):
    """Calculate the number of days between a given date and a reference day-month (year agnostic, cyclic)
    Args:
        date (datetime): single value or list of datetimes, for calculating against reference date. 
        ref_month (str): reference month
        ref_day (str): reference day
    Returns:
        integer: days since reference day-month (value will be between 0-365/366).
    """
    date = pd.to_datetime(date)
    ref_date = pd.to_datetime(f'{date.year}-{ref_month}-{ref_day}',format='%Y-%m-%d')
    
    if date < ref_date:
        ref_date = ref_date - pd.DateOffset(years=1)
    
    return (date - ref_date).days

def cumsum_reset_at_value(da,reset_value=0):
    """Apply a cumulative sum along the time axis, with cumsum reset when encounter reset_value
    Args:
        da (xarray dataarray): dataarray to apply cumsum to
        reset_value (float): when this value is encountered, cumsum reset to zero
    Returns:
        xarray dataarray: array with cumulative summed values (reset as appropriate)
    """    
    without_reset = da.values.cumsum(axis=da.get_axis_num('time'))
    reset_at = (da.values == reset_value)
    overcount = np.maximum.accumulate(without_reset * reset_at)
    cumsum = without_reset - overcount
    return xr.DataArray(data=cumsum,dims=da.dims,coords=da.coords)

def convert_calendar(ds, output_calendar):
    """Convert time calendar.
    Original Author: Damien Irving
    """

    valid_calendars = {
        cftime._cftime.DatetimeGregorian: cftime.DatetimeGregorian,
        cftime._cftime.DatetimeProlepticGregorian: cftime.DatetimeProlepticGregorian,
        cftime._cftime.DatetimeNoLeap: cftime._cftime.DatetimeNoLeap,
    }

    output_calendar_name = str(output_calendar).split('.')[-1][:-2]
    if output_calendar in valid_calendars:
        input_calendar_name = str(type(ds['time'].values[0])).split('.')[-1][:-2]
        output_calendar_name = str(output_calendar).split('.')[-1][:-2]
        logging.info(f'Convering input {input_calendar_name} calendar to {output_calendar_name}')

        is_noleap = output_calendar == cftime._cftime.DatetimeNoLeap
        if is_noleap:
            ds = ds.sel(time=~((ds['time'].dt.month == 2) & (ds['time'].dt.day == 29)))

        new_times = []
        calendar_func = valid_calendars[output_calendar]
        for old_time in ds['time'].values:
            new_time = calendar_func(old_time.year, old_time.month, old_time.day, old_time.hour)
            new_times.append(new_time)
        time_attrs = ds['time'].attrs
        ds = ds.assign_coords({'time': new_times})
        ds['time'].attrs = time_attrs

        if 'time_bnds' in ds:
            new_time_bnds = []
            for old_start, old_end in ds['time_bnds'].values:
                if is_noleap and (old_start.day == 29) and (old_start.month == 2):
                    old_start_day = 28
                    old_start_month = 2
                else:
                    old_start_day = old_start.day
                    old_start_month = old_start.month
                if is_noleap and (old_end.day == 29) and (old_end.month == 2):
                    old_end_day = 1
                    old_end_month = 3
                else:
                    old_end_day = old_end.day
                    old_end_month = old_end.month
                new_start = calendar_func(old_start.year, old_start_month, old_start_day, old_start.hour)
                new_end = calendar_func(old_end.year, old_end_month, old_end_day, old_end.hour)
                time_diff = new_end - new_start
                assert time_diff == np.timedelta64(1, 'D')
                new_time_bnds.append([new_start, new_end])

            da_time_bnds = xr.DataArray(
                new_time_bnds,
                dims=ds['time_bnds'].dims,
                coords={"time": ds['time']},
            )
            ds['time_bnds'] = da_time_bnds
    else:
        raise ValueError(f'Conversion to {output_calendar_name} not supported')

    return ds

# The following functions are related to adding time_bounds to datasets

def _get_new_time_bnds(time_bnds, resamp):
    """Get time bounds using old time bounds and the resampling.
    Args:
        time_bnds (list or xarray): time bounds
        resamp (dict): resampling (e.g., {time:1D})
    Returns:
        xarray : time bounds
    """
    new_time_bnds = []
    time_bnds = time_bnds.load() # Makes it faster
    for i in time_bnds.resample(resamp)._iter_grouped():
        new_time_bnds.append( [i.isel(bnds=0, time=0).values, i.isel(bnds=-1, time=-1).values] ) #< Assuming time is sorted select the first and last time bounds in each resample group as new bounds
    return new_time_bnds

def get_new_time_bnds(ds, ds_out, resamp):
    """Get time bounds using old new and old dataset and the resampling.
    Args:
        ds (xarray dataset): Original dataset before resampling
        ds_out (xarray dataset): New dataset after resampling
        resamp (dict): resampling (e.g., {time:1D})
    Returns:
        xarray : time bounds
    """
    if "time_bnds" in ds:
        new_time_bnds       = _get_new_time_bnds(ds["time_bnds"], resamp)
        ds_out['time_bnds'] = (("time", "bnds"), new_time_bnds)
    return ds_out

def guess_time_bnds(ds, ds_out, resamp):
    """Guess time bounds using old new and old dataset and the resampling.
    Args:
        ds (xarray dataset): Original dataset before resampling
        ds_out (xarray dataset): New dataset after resampling
        resamp (dict): resampling (e.g., {time:1D})
    Returns:
        xarray : time bounds
    """
    dt_unit = list(resamp.values())[0]

    try:
        lower_bound = ds['time'].resample(resamp).min().dt.floor(dt_unit)
    except ValueError:
        #< dt.floor does not work with non-fixed frequency (e.g., month). Thus floor to the closest day instead
        dt_unit = "D"
        lower_bound = ds['time'].resample(resamp).min().dt.floor(dt_unit)
        
    upper_bound = ds['time'].resample(resamp).max().dt.ceil(dt_unit)


    new_time_bnds = []
    for i in range(len(ds_out['time'])):
        new_time_bnds.append( [lower_bound[i].values, upper_bound[i].values] )

    ds_out['time_bnds'] = ds_out['time'].expand_dims(dim={'bnds':2})
    ds_out['time_bnds'] = (("time", "bnds"), new_time_bnds)
    return ds_out

def update_time_bnds(ds, ds_out, resamp):
    """Update time bounds using old new and old dataset and the resampling.
    Args:
        ds (xarray dataset): Original dataset before resampling
        ds_out (xarray dataset): New dataset after resampling
        resamp (dict): resampling (e.g., {time:1D})
    Returns:
        xarray : time bounds
    """
    if "time_bnds" in ds:
        ds_out = get_new_time_bnds(ds, ds_out, resamp)
    else:
        ds_out = guess_time_bnds(ds, ds_out, resamp)

    #< Place the valid time in the middle of the time bounds
    ds_out["time"] = ds.time.resample(resamp).mean()

    # Fix an issue where cftime bounds are arrays of cftime instead of cftime
    if isinstance(ds_out["time_bnds"].values[0][0], np.ndarray):
        for i in range(len(ds_out["time_bnds"])):
            ds_out["time_bnds"][i][0] = ds_out["time_bnds"][i][0].values[()]
            ds_out["time_bnds"][i][1] = ds_out["time_bnds"][i][1].values[()]
    
    return ds_out

def center_valid_time_between_time_bnds(ds):
    """If time bounds are present in the dataset make sure the valid time is between the time_bnds
    Args:
        ds (xarray dataset): Dataarray to check valid time is centred
    Returns:
        xarray dataset : New dataset with centred valid time
    """
    #< If time_bnds is not present exit
    if not "time_bnds" in ds:
        return ds

    else:
        ds["time"] = ds["time_bnds"][:,0] + (ds["time_bnds"][:,1] - ds["time_bnds"][:,0])/2
    return ds
 
def add_attrs(arr, key, value):
    """Add attributes to a xarray dataarray. If attribute exists add to existing attribute.
    Args:
        arr (xarray dataarray): Dataarray to add attribute to
        key (object): Key of the new attribute
        value (object): Value of the new attribute
    Returns:
        xarray : New array with attribute added
    """
    if key in arr.attrs:
        arr.attrs[key] += " " + value
    else:
        arr.attrs[key] = value
    return arr

def _add_cell_methods(arr, cell_method):
    """Add cell method to a xarray dataarray. If cell method exists add to existing.
    Args:
        arr (xarray dataarray): Dataarray to add cell methods to
        cell_method (str): cell_method (e.g. time: mean)
    Returns:
        xarray : New array with cell_method added
    """
    return add_attrs(arr, key="cell_methods", value=cell_method)

def add_cell_methods(ds, resamp=None, dim=None, operation="", exclude_list=[]):
    """Add cell method to a xarray dataset being resampled. If cell method exists add to existing.
    Args:
        ds (xarray dataset): Dataset to add cell methods to
        resamp (dict): resampling (e.g., {time:1D})
        dim (str): If not resamp is provided but dim add cell method for operation over complete dimension.
        operation (str): Operation to apply (e.g., mean max min)
        exclude_list (list): Variables in this list get excluded
    Returns:
        xarray : New dataset with cell_method added
    """
    for v in ds: # Go through all variables in dataset
        if not v in exclude_list: #< Only if not excluded
            if resamp is None: #< If there is not resampling check if a dimension was defined
                if not dim is None: #< If dim is defined use it in cell_methods
                    cell_method = f"{dim}: {operation}"
                    ds[v]       = _add_cell_methods(ds[v], cell_method)
                else:   #< Else the operation is assumed to go over all dimensions
                    dims        = ": ".join(ds.dims)
                    cell_method = f"{dims}: {operation}"
                    ds[v]       = _add_cell_methods(ds[v], cell_method)
            else:
                for r in resamp:
                    cell_method = f"{r}: {operation} (interval: {resamp[r]})" if resamp is not None else f"{r}: {operation}"
                    ds[v]       = _add_cell_methods(ds[v], cell_method)
    return ds



