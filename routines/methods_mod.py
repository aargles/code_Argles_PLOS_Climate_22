# -*- coding: utf-8 -*-
"""
Created on 22/04/2022

Defined here are the functions used in making the data for figures 1 and 2 in Argles, et al. in review. Output data is saved in the ‘data’ subfolder in the principal folder. Specifically, estimation of the minimum and maximum grid-box areas for a given longitude and latitude grid resolution (Fig 1), and simulation of two forest patch-models (Fig 2). 

The script utilises definitions defined locally in “methods_submod.py”. In there is defined the actual numerical simulation of the patch dynamics.

References
----------

Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

"""

# Credits
__author__ = "Arthur Argles, Jonathan Moore, and Peter Cox"
__credits__ = ["Arthur Argles","Jonathan Moore","Peter Cox"]
__license__ = "CCBY4.0"
__maintainer__ = "Arthur Argles"
__email__ = "A.Argles2@exeter.ac.uk"

import pickle
import numpy as np
import pandas as pd


from .methods_submod import *


def generate_fig_1_gridbox_area(data_dir,datafile='fig01.csv',lon_res=1.88,\
                                lat_res=1.25):
    """
    Generates and saves data used in fig 1 in Argles, et al. (in review). This function calculates the area of a grid-box for fig01.
  
    References
    ----------
    
    Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

    Parameters
    ----------
    data_dir : str
        Where the model output data is stored.
    datafile : str, optional
        The filename of the output data file.
    lon_res : float, optional
        The longitudinal resolution of the grid. The default is 1.88. (degrees)
    lat_res : TYPE, optional
        The latitudinal resolution of the grid. The default is 1.25. (degrees)

    Returns
    -------
    None.

    """
    
    ### Data
    # datafile Format
    # ---------------
    #  For the correct formatting, datafile should be a csv file with the column headers: "feature name", "space center (m2)", "time center (yr)", "min space (m2)", "max space (m2)", "min time (yr)", "max time (yr)", "plot feature", "text pos", "text rotation", "zorder"
    #
    #  Column are defined in data_dir/readme.txt under the fig01 section
    #
    #  For this function to properly work we need to have "Gridbox Area Min", and "Gridbox Area Max" already defined, as we are just modifying the file
    veg_scale_df = pd.read_csv(data_dir+'/'+datafile)
    
    ### Estimate area
    # The default value is assumed to follow the N96 grid of the UK Met Office Unified Model Standard
    area_min, area_max = get_max_min_area_on_earth_grid(lon_res=lon_res,lat_res=lat_res)
    
    i_gridbox_min = \
        np.squeeze(np.argwhere(np.array(veg_scale_df['feature name']==\
                                        'Gridbox Area Min')))
    i_gridbox_max = \
        np.squeeze(np.argwhere(np.array(veg_scale_df['feature name']==\
                                        'Gridbox Area Max')))
    # Overwrite saved file with the estimated grid-box area
    veg_scale_df.iloc[i_gridbox_min].at['min space (m2)'] = area_min
    veg_scale_df.iloc[i_gridbox_min].at['max space (m2)'] = area_min
    veg_scale_df.iloc[i_gridbox_max].at['min space (m2)'] = area_max
    veg_scale_df.iloc[i_gridbox_max].at['max space (m2)'] = area_max
    veg_scale_df.to_csv(data_dir+'/'+datafile,index=False)
    
    return
    
    
    

def generate_fig_2_patch_data(t_end,dt,loss_rate,data_dir,\
                              datafile='fig02.pickle',include_VF=True,\
                              include_grid=True,grid_res=[50,50],\
                              t_slices=np.array([5.,50.,250.,500.])):
    """
    Generates and saves data used in fig 2 in Argles, et al. (in review). This function simulates two possible methods of capturing the horizontal distribution of forests in the form of patch age assuming a constant rate of disturbance.
  
    References
    ----------
    
    Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).
    
    Parameters
    ----------
    t_end : float
        The total number of years for the simulation (with t_start = 0). (yr)
    dt : float
        The timestep of the simulation. (yr)
    loss_rate : float
        The rate of patch disturbance. (/yr)
    data_dir : str
        Where the model output data is stored.
    datafile : str, optional
        The filename of the output data file. The default is 'fig02.pickle'.
    include_VF : bool, optional
        If true simulate the Von Foerster patch-dynamics equation via the methods of characteristics. The default is True.
    include_grid : bool, optional
        If true simulate a random grid of patches of resolution grid_res. The default is True.
    grid_res :  array_like[size=(2,),dtype=int], optional
        The grid-resolution of the random run, the total number of patches simulated is the product of the elements. The default is [50,50].
    t_slices :  array_like[size=(N,),dtype=float], optional
        Where we want to save the spatial and age-structure time-slices for fig 2. The default is np.array([5.,50.,250.,500.])
        
    Raises
    ------
    UserWarning
        If no model boolean is set to True then we cannot run and save any data.

    Returns
    -------
    None.

    """
        
    # Check to see if we are running at least one model.    
    if (include_VF == False and include_grid == False):
    
        raise UserWarning('No patch data will be simulated. For running '+\
                          'the toy models please set at least one model '+\
                          'boolean to True (e.g. inclue_VF = True).')
        
    # Get the number of slices for age-structure as output
    K_slice = len(t_slices)
    
    # Get duration of simulation and the number of timesteps
    t_start = 0.
    t = np.append(np.arange(t_start,t_end,step=dt),t_end)
    K = len(t)

    # Find closest elements next to time slices
    k_sclices = np.zeros(K_slice,dtype=int)
    
    for k in range(0,K_slice):
        
        k_sclices[k] = find_nearest(t,t_slices[k])

    # For saving into output pickle file as we go along
    output_dict = {'t_slices':t_slices,'k_sclices':k_sclices,\
                   'grid_res':grid_res,'t':t,'dt':dt,'loss_rate':loss_rate}
        
    ### Von Foerster patch-dynamics
    if include_VF == True:
        
        # Like some DGVMs such as ED (Moorecraft, et al., 2001) we use the methods of characteristics
        patch_age_pdf = patch_VF_MOC_dynamic(K,dt,loss_rate)
        # assume dt = dage
        age_VF = np.zeros(K)
        age_VF[:] = t[:]
        # For the number of patches with positive pdf over time 
        num_patches_VF = np.zeros(K)
        num_patches_VF = \
            np.array([len(np.argwhere(patch_age_pdf[k,:]>0)) for k in\
                      range(0,K)],dtype=int)
        
        patch_age_pdf_slice = np.zeros(K_slice,dtype=object)
        
        # Run through the time slices to record the patch distribution wrt patch age
        for k in range(0,K_slice):
            
            patch_age_pdf_slice[k] = patch_age_pdf[k_sclices[k]]
                                    
        # Here we update our output dict to save the data
        output_dict['age_VF'] = age_VF
        output_dict['patch_age_pdf_VF'] = patch_age_pdf_slice
        output_dict['num_patches_VF'] = num_patches_VF
        
    ### Stochastic patch-dynamics
    if include_grid == True:
                
        I = grid_res[0]
        J = grid_res[1]
        patch_age_stochastic = patch_grid_stochastic(K,I,J,dt,loss_rate)
        patch_age_stochastic_grid = np.zeros((K_slice,I,J))
        
        # We allow the bins and number of bins to vary across the number of slices, this is so we can clearly demonstrate the change in distribution,
        age_bins_S = np.zeros(K_slice,dtype=object)
        age_bin_centers_S  = np.zeros(K_slice,dtype=object)
        age_bin_widths_S = np.zeros(K_slice,dtype=object)
        patch_age_pdf_S = np.zeros(K_slice,dtype=object)
        
        # For loop runs through each of the time-slices recording the spatial patch-layout and the patch distribution wrt patch age
        for k in range(0,K_slice):
            
            patch_age_stochastic_grid[k,:,:] = \
                patch_age_stochastic[k_sclices[k],:,:]
            patch_age_stochastic_sclices_flat = \
                np.ndarray.flatten(patch_age_stochastic_grid[k,:,:])
            age_bins_S[k], \
            age_bin_centers_S[k], \
            age_bin_widths_S[k] = \
                equiprobable_bins(patch_age_stochastic_sclices_flat)         
            patch_age_pdf_S[k], _ =\
                np.histogram(patch_age_stochastic_sclices_flat,\
                bins=age_bins_S[k],density=True)
                
        # Here we update our output dict to save the data
        output_dict['patch_age_stochastic_grid'] = patch_age_stochastic_grid
        output_dict['age_bins_S'] = age_bins_S
        output_dict['age_bin_centers_S'] = age_bin_centers_S
        output_dict['age_bin_widths_S'] = age_bin_widths_S
        output_dict['patch_age_pdf_S'] = patch_age_pdf_S        
    
    ### Save Data
    with open(data_dir+'/'+datafile, 'wb') as handle:
        pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    return
    
    
    