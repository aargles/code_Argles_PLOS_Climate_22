# -*- coding: utf-8 -*-
"""
Created on 22/04/2022

Defined here are definitions that are used locally in "methods_mod.py”, these definitions are primarily to keep code tidy and save space.

Internally, under section "#%% Primary Definitions" we have the simulations and calculation of the stated objectives: "get_max_min_area_on_earth_grid" and simulation of patch dynamics statistically and stochastically.

Under section "#%% Secondary Definitions" we have defined practical functions (e.g. find the nearest vector index to a target value) to save on space. 

"""

# Credits
__author__ = "Arthur Argles, Jonathan Moore, and Peter Cox"
__credits__ = ["Arthur Argles","Jonathan Moore","Peter Cox"]
__license__ = "CCBY4.0"
__maintainer__ = "Arthur Argles"
__email__ = "A.Argles2@exeter.ac.uk"


import numpy as np
from math import sin, cos, sqrt, atan2, radians


#%% Primary Definitions

def get_max_min_area_on_earth_grid(lon_res=1.88,lat_res=1.25):
    """
    Function estimates the maximum (near equator) and minimum (near poles) grid-box area (m2) based around the inputted longitude and latitude resolutions. This method assumes a spherical earth.

    Parameters
    ----------
    lon_res : float, optional
        The longitudinal resolution of the grid. The default is 1.88. (degrees)
    lat_res : TYPE, optional
        The latitudinal resolution of the grid. The default is 1.25. (degrees)

    Returns
    -------
    area_min : float
        The minimum grid-box area of the grid. Assumed to be closest grid-box to 90 degrees north. (m2)
    area_max : float
        The maximum grid-box area of the grid. Assumed to be closest grid-box to 0 degrees north. (m2)

    """
    
    # Taking -180 to 180 degree west for longitude limits and -90 degree to 90 degree north for latitude limits.
    lon = np.arange(-180+lon_res/2.0,180-lon_res/2.0,step=lon_res)
    lat = np.arange(-90+lat_res/2.0,90-lat_res/2.0,step=lat_res)
    n_lon = len(lon)
    m_lat = len(lat)
    
    # Make corresponding grid
    lon_xx, lat_yy = np.meshgrid(lon,lat)
    
    # Run through each gridpoint and estimate the grid-box area, add up to check total area is accurate (~5.1e14 m2)
    area = np.zeros((m_lat,n_lon))
    total_area = 0.
    
    for n in range(0,n_lon):
        
        for m in range(0,m_lat):
            
            area[m,n] = quadrant_lonlat_area(lon[n],lat[m],lon_res,lat_res)
            total_area = area[m,n] + total_area    
    
    # For checking method       
    #print("{:e}".format(total_area))
            
    area = area[:,0]
    degree_90_elm = find_nearest(lat,90)
    degree_0_elm = find_nearest(lat,0)
    area_min = area[degree_90_elm]
    area_max = area[degree_0_elm]
    
    return area_min, area_max

def patch_VF_MOC_dynamic(K,dt,loss_rate,patch_age_pdf_0=None):
    """
    Simulate the patch-age Probablity Density Function using the Methods of Characteristics.    

    Parameters
    ----------
    K : int
        The total number of time-steps.
    dt : float
        The timestep of the simulation. (yr)
    loss_rate : float
        The rate of patch disturbance. (/yr)
    patch_age_pdf_0 : array_like[size=(K,),dtype=float], optional
        Thye initial patch age pdf distribution. The default (patch_age_pdf_0=None) assumes everything starts at age=0.

    Returns
    -------
    patch_age_pdf : array_like[size=(K,K),dtype=float]
        The evolution of the patch age pdf through simulation time. (/yr)

    """
         
    if patch_age_pdf_0 is None:
        
        patch_age_pdf_0  = np.zeros(K)
        patch_age_pdf_0[0] = 1.0/dt  # As we assume that da = dt for the Methods of Characteristics
        
    patch_age_pdf = np.zeros((K,K))
    patch_age_pdf[0,:] = patch_age_pdf_0[:]
    
    for k in range(1,K):
        
        for k_age in range(0,K):
            
            if k_age > 0:
                
                patch_age_pdf[k,k_age] = \
                    patch_age_pdf[k-1,k_age-1] * np.exp(-loss_rate*dt)
                
            else:
                
                patch_age_pdf[k,0] = \
                    np.sum(patch_age_pdf[k-1,:]*(1.0 - np.exp(-loss_rate*dt)))
                    
    return patch_age_pdf

def patch_grid_stochastic(K,I,J,dt,loss_rate,patch_age_0=None):
    """
    Run a stochastic toy model of patch-age assuming a constant rate of patch disturbance.

    Parameters
    ----------
    K : int
        The total number of time-steps.
    I : int
        The total number of patches along the i axis.
    J : int
        The total number of patches along the j axis.
    dt : float
        The timestep of the simulation. (yr)
    loss_rate : float
        The rate of patch disturbance. (/yr)
    patch_age_0 : array_like, optional
        The initial patch age distribution of the grid. The default is None.

    Returns
    -------
    patch_age : array_like[size=(K,I,J),dtype=float]
        The evolution of the gridded patches through time.(/yr)

    """
    
    # Just in case we want to start with an initial distribution
    if patch_age_0 is None:
        
        patch_age_0 = np.zeros((I,J))
        
    patch_age = np.zeros((K,I,J))
    patch_age[0,:,:] = patch_age_0[:,:]
    
    
    # The chance of death is an exponential curve across the time-step, when integrated. This is because of the constant loss_rate.
    chance = 1.0 - np.exp(-loss_rate*dt)
    nochance = 1.0 - chance
    # To speed things up we use the random.choice function to select when a patch dies across space and time.
    mort = np.random.choice([1,0],size=(K-1,I,J),p=[chance,nochance])
    
    # For loop runs through time adding each timestep, or if patch has died reseting patch-age to zero
    for k in range(1,K):
        
        patch_age[k,:,:] = patch_age[k-1,:,:] + dt
        patch_age[k,mort[k-1,]==1] = 0.
    
    return patch_age


#%% Secondary Definitions

def find_nearest(array, value):
    """
    From a target value find the index of the nearest value from an array.    
    
    Parameters
    ----------
    array : array_like[size=(N,),dtype=float]
        Array we are finding index of nearest value to target.
    value : float
        Target value.

    Returns
    -------
    idx : int
        Index of nearest value.

    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def lonlat_to_m(point_1,point_2):
    """
    Spherical earth distance in meters between two longitude and latitude points

    Parameters
    ----------
    point_1 : array_like[size=(2,),dtype=float]
        Coordinantes of the first point in longitude and latitude (lon,lat). (degrees)
    point_2 : array_like[size=(2,),dtype=float]
        Coordinantes of the second point in longitude and latitude (lon,lat). (degrees)

    Returns
    -------
    d : float
        Arc distance between points, assuming the earth is a spehere. (m)

    """
    
    R = 6.378e6
    lon_1 = radians(point_1[0])
    lon_2 = radians(point_2[0])
    lat_1 = radians(point_1[1])
    lat_2 = radians(point_2[1])
    dlon = lon_2 - lon_1
    dlat = lat_2 - lat_1
    a = sin(dlat/2)**2+cos(lat_1)*cos(lat_2)*sin(dlon/2)**2
    c = 2*atan2(sqrt(a),sqrt(1-a))  
    d = R*c
    
    return d
    

def quadrant_lonlat_area(lon,lat,lon_res,lat_res):
    """
    For a grid-box on assumed spherical Earths Surface estimate the area in m2 based on longitude, latitude, and resolution. 

    Parameters
    ----------
    lon : float
        The longitude of the grid-box. (degrees)
    lat : float
        The latitude of the grid-box (degrees).
    lon_res : float
        The longitudinal resolution of the grid-box. (degrees)
    lat_res : float
        The latitudinal resolution of the grid-box. (degrees)

    Returns
    -------
    area : float
        The area of the grid-box. (m2)

    """
    
    # Get the boundaries of the grid-box (assumed to be halfway)
    point_lower = (lon,lat-lat_res/2.0)
    point_upper = (lon,lat+lat_res/2.0)
    point_left = (lon+lon_res/2.0,lat)
    point_right = (lon-lon_res/2.0,lat)
    
    # Get the height and width of the grid-box in m2, then estimate the area
    height = lonlat_to_m(point_lower,point_upper)
    width = lonlat_to_m(point_left,point_right)
    area = height * width
    
    return area



def equiprobable_bins(data):
    """        
    Function generates equiprobable bins, bin_centers and bin_widths, for the input data.


    Parameters
    ----------
    data : array_like[size=(N,),dtype=float]
        The input data we are proscribing bins to.

    Returns
    -------
    bins : array_like[size=(K,),dtype=float]
        The equiprobable bins for the given data.
    bin_centers : array_like[size=(K-1,),dtype=float]
        The equiprobable bin centers for the given data.
    bin_widths : array_like[size=(K-1,),dtype=float]
        The equiprobable bin widths for the given data.
        
    """
    
    # Broad rule based on https://itl.nist.gov/div898/handbook/prc/section2/prc211.htm
    K = int(round(2 * (len(data))**(2/5)))
    percentile_bins = np.linspace(0,100,K)
     
    bins = np.zeros(K)
    
    # For loop runs through each percentile.
    for k in range(0,K):
        
        bins[k] = np.percentile(data,percentile_bins[k])
        
    # Eliminate repeat bin values (only necessary if data is near homogeneous)
    bins = np.unique(bins)
    K = len(bins)
    bin_centers = np.zeros(K-1)
    bin_widths = np.zeros(K-1)
    
    for k in range(1,K):
        
        bin_widths[k-1] = bins[k] - bins[k-1]
        bin_centers[k-1] = bins[k-1] + 0.5*bin_widths[k-1]
    
    return bins, bin_centers, bin_widths