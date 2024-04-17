""" Functions related to reading and writing files """

import xarray as xr
import pandas as pd
import numpy as np  
from datetime import date
from datetime import datetime
import calendar
import sys, os
import glob
import git
import cmdline_provenance as cmdprov
import xesmf as xe
import xclim as xc
import logging

def read_data(
    infiles,
    var,
    time_bounds=None,
    lat_bounds=None,
    lon_bounds=None,
    input_units=None,
    output_units=None,
    lon_chunk_size=None,
    apply_ssr=False,
    use_cftime=True,
    output_calendar=None,
    valid_min=None,
    valid_max=None,
    return_as_DataArray=True,
):
    """Read and process an input dataset.

    Parameters
    ----------
    infiles : list
        Input files    
    var : str
        Variable to read
    time_bounds : list, optional
        Time period to extract from infiles [YYYY-MM-DD, YYYY-MM-DD]
    lat_bnds : list, optional
        Latitude bounds: [south bound, north bound] 
    lon_bnds : list, optional
        Longitude bounds: [west bound, east bound]    
    input_units : str, optional
        Units of input data (if not provided will attempt to read file metadata)
    output_units : str, optional
        Desired units for output data (conversion will be applied if necessary)
    lon_chunk_size : int, optional
        Put this number of longitudes in each data chunk
    apply_ssr : bool, default False
        Apply Singularity Stochastic Removal to the data
    use_cftime : bool, default True
        Use cftime for time axis
    output_calendar : cftime calendar, optional
        Desired calendar for output data
    valid_min : float, optional
        Clip data to valid minimum value
    valid_max : float, optional
        Clip data to valid maximum value
    return_as_DataArray : bool, default True
        Returns xarray DataArray instead of DataSet

    Returns
    -------
    ds : xarray DataArray / DataSet

    Original author: Damien Irving

    """

    if len(infiles) == 1:
        try:
            ds = xr.open_dataset(infiles[0], use_cftime=use_cftime)
        except ValueError:
            ds = xr.open_dataset(infiles[0])
    else:
        try:
            ds = xr.open_mfdataset(infiles, use_cftime=use_cftime)
        except ValueError:
            ds = xr.open_mfdataset(infiles)

    try:
        ds = ds.drop('height')
    except ValueError:
        pass

    if 'latitude' in ds.dims:
        ds = ds.rename({'latitude': 'lat'})
    if 'longitude' in ds.dims:
        ds = ds.rename({'longitude': 'lon'})

    if time_bounds:
        start_date, end_date = time_bounds
        ds = ds.sel({'time': slice(start_date, end_date)})        
    if lat_bounds:
        ds = subset_lat(ds, lat_bounds)
    if lon_bounds:
        ds = subset_lon(ds, lon_bounds)

    if output_calendar:
        input_calendar = type(ds['time'].values[0])
        if input_calendar != output_calendar:
            ds = convert_calendar(ds, output_calendar)  

    if input_units:
        ds[var].attrs['units'] = input_units
    if output_units:
        ds[var] = convert_units(ds[var], output_units)
        ds[var].attrs['units'] = output_units

    if (valid_min is not None) or (valid_max is not None):
        ds[var] = ds[var].clip(min=valid_min, max=valid_max, keep_attrs=True)

    chunk_dict = {'time': -1}
    if lon_chunk_size:
        chunk_dict['lon'] = lon_chunk_size
    ds = ds.chunk(chunk_dict)
    
    if return_as_DataArray:
        ds = ds[var]

    return ds

def save_data(ds,ofile,provenance=True):
    """Write an xarray DataSet/DataArray to file

    Parameters
    ----------
    ds : xarray DataSet or DataArray
        Data to write to file
    ofile : str
        Output file name
    provenance : str
        Add commandline provenance
    """
    
    if provenance:
        ds.attrs['history'] = cmdprov.new_log(extra_notes=[get_git_hash()])
    
    ds.to_netcdf(ofile)

def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text

def assert_is_DataArray(da):
    """Assert input data is of type xarray DataArray"""

    assert isinstance(da, xr.core.dataarray.DataArray), f"input data needs to be an xarray DataArray, is type {type(da)}"

def assert_is_Dataset(ds):
    """Assert input data is of type xarray Dataset"""

    assert isinstance(ds, xr.core.dataset.Dataset), f"input data needs to be an xarray Dataset, is type {type(ds)}. Consider using da.to_dataset()"

def update_attrs(da, attrs):
    """Update attributes of xarray DataArray

    Parameters
    ----------
    da : xarray DataArray
        Input array to update attributes
    attrs : dict
        Dictionary containing attributes to update.

    Returns
    -------
    da : xarray DataArray
        Array with updated attributes
    """

    assert_is_DataArray(da)
    
    if 'name' in attrs.keys():
        da.name= attrs['name']
        del attrs['name']
    for a in attrs.keys():
        da.attrs[a] = attrs[a]
    
    return da

def convert_units(da, target_units):
    """Convert units.

    Parameters
    ----------
    da : xarray DataArray
        Input array containing a units attribute
    target_units : str
        Units to convert to

    Returns
    -------
    da : xarray DataArray
       Array with converted units
    
    Original author: Damien Irving
    """

    custom_conversions = {
        ("MJ m-2", "W m-2"): joules_to_watts,
        ("megajoule/meter2", "W m-2"): joules_to_watts,
    }
    try:
        da = xc.units.convert_units_to(da, target_units)
    except Exception as e:
        var_attrs = da.attrs
        conversion = (da.attrs["units"], target_units)
        if conversion in custom_conversions:
            da = custom_conversions[conversion](da)
            da.attrs = var_attrs
        else:
            raise e

    return da

def joules_to_watts(da):
    """Convert from Joules to Watts"""

    input_units = da.attrs["units"]
    input_freq = xr.infer_freq(da.indexes['time'][0:3])[0]
    assert input_freq == 'D'

    if (input_units[0] == 'M') or (input_units[0:4] == 'mega'):
        da = da * 1e6
    seconds_in_day = 60 * 60 * 24
    da = da / seconds_in_day

    return da

def apply_ssr(da, threshold='8.64e-4 mm day-1'):
    """Apply Singularity Stochastic Removal.

    Used to avoid divide by zero errors in the analysis of precipitation data.
    All near-zero values (i.e. < threshold) are set to a small random non-zero value:
    0 < value <= threshold
    
    Parameters
    ----------
    da : xarray DataArray
        Input precipitation data
    threhsold : str, default '8.64e-4 mm day-1'
        Threshold for near-zero rainfall

    Returns
    -------
    da_ssr : xarray DataArray
        Input data with ssr applied

    Reference
    ---------
    Vrac, M., Noel, T., & Vautard, R. (2016). Bias correction of precipitation
    through Singularity Stochastic Removal: Because occurrences matter.
    Journal of Geophysical Research: Atmospheres, 121(10), 52375258.
    https://doi.org/10.1002/2015JD024511

    Original Author: Damien Irving
    """

    da_ssr = sdba.processing.jitter_under_thresh(da, '8.64e-4 mm day-1')

    return da_ssr


def reverse_ssr(da_ssr, threshold=8.64e-4):
    """Reverse Singularity Stochastic Removal.

    SSR is used to avoid divide by zero errors in the analysis of precipitation data.
    It involves setting near-zero values (i.e. < threshold) to a small non-zero random value: 0 < value <= threshold.
    This function reverses SSR (commonly at the end of a calculation) by setting all near-zero values (i.e. < threshold) to zero.
    
    Parameters
    ----------
    da_ssr : xarray DataArray
        Input precipitation data (that has had SSR applied)
    threhsold : float, default 8.64e-4 mm
        Threshold for near-zero rainfall

    Returns
    -------
    da_no_ssr : xarray DataArray
        Input data with ssr reversed

    Reference
    ---------
    Vrac, M., Noel, T., & Vautard, R. (2016). Bias correction of precipitation
    through Singularity Stochastic Removal: Because occurrences matter.
    Journal of Geophysical Research: Atmospheres, 121(10), 52375258.
    https://doi.org/10.1002/2015JD024511
    
    Original Author: Damien Irving
    """

    da_no_ssr = da_ssr.where(da_ssr >= threshold, 0.0)

    return da_no_ssr
