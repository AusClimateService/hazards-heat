import numpy as np
import xarray as xr
import pandas as pd

def get_GWL_syear_eyear(CMIP,GCM,pathway,GWL):
    """Returns the start and end year of the Global Warming Level timeslice for the specified GWL, GCM and pathway.
    This script searches the IPCC-WG1 Atlas lookup table of Global Warming Levels:\
    https://github.com/IPCC-WG1/Atlas/tree/main/warming-levels

    Parameters
    ----------
    CMIP : str
        Version of CMIP [options: CMIP5, CMIP6]
    GCM : str
        Name of the global climate model
    pathway : str
        Emissions pathway that GCM run under[options: rcp26, rcp45, rcp85, ssp126, ssp245, ssp370, ssp585]
    GWL : str
        Required Global Warming Level [options: current, 1.5, 2, 3, 4]

    Returns
    -------
    int : syear, eyear
       Start and end year of corresponding 20-year GWL timeslice
    """

    pathway = pathway.lower()
    assert pathway in ['rcp26','rcp45','rcp85','ssp126','ssp245','ssp370','ssp585']
    assert GWL in ['current','1.5','2','3','4']

    if GWL == 'current':
        gwl_center = 2020
    else:
        df = pd.read_csv('https://raw.githubusercontent.com/IPCC-WG1/Atlas/main/warming-levels/{}_Atlas_WarmingLevels.csv'.format(CMIP))
        df = df.replace(9999, np.NaN)
        df['model_run'] = df['model_run'].str.lower()
        gwl_center = df[df['model_run'].str.contains(GCM.lower())][GWL+"_"+pathway].values[0]
    syear = int(gwl_center)-9
    eyear = int(gwl_center)+10

    print(GCM, GWL, gwl_center, syear,'-', eyear)
    return syear, eyear

def get_GWL_timeslice(ds,CMIP,GCM,pathway,GWL):
    """Returns the 20-year timeslice of a data array corresponding to the desired Global Warming Level.

    Parameters
    ----------
    ds : xarray DataArray or DataSet
        Input array from which to take the 20-year timeslice
    CMIP : str
        Version of CMIP [options: CMIP5, CMIP6]
    GCM : str
        Name of the global climate model
    pathway : str
        Emissions pathway that GCM run under[options: rcp26, rcp45, rcp85, ssp126, ssp245, ssp370, ssp585]
    GWL : str
        Required Global Warming Level [options: current, 1.5, 2, 3, 4]

    Returns
    -------
    ds : xarray DataArray or DataSet
       Subsampled data array corresponding to the desired Global Warming Level.
    """

    syear,eyear = get_GWL_syear_eyear(CMIP,GCM,pathway,GWL)
    return ds.sel(time=slice('{}-01-01'.format(int(syear)),'{}-12-31'.format(int(eyear))))

def emission_pathway(Y,pathway):
    """Reset pathway to 'historical' for years earlier than start year of SSP/RCP runs
    Args:
        Y (int): year
        pathway (str): future emission pathway
    Returns:
        string: appropriate model pathway for given model year.
    """
    if pathway.lower().replace('.','') in ['ssp126','ssp245','ssp370','ssp585']:
        if Y < 2015:
            pathway = 'historical'
    elif pathway.lower().replace('.','') in ['rcp26','rcp45','rcp85']:
        if Y < 2005:
            pathway = 'historical'
    elif pathway != 'historical':
        raise ValueError(f"Specified pathway {pathway} not recognised")
    
    return pathway


