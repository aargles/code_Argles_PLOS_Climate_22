# -*- coding: utf-8 -*-
"""
Created on 22/04/2022

Defined here are the functions used in making figures 1, 2 and 3 in Argles et al., in review. The figure captions are declared in the function descriptions. Data for these figures is stored in the ‘data’ subfolder below the principal folder. 

The python script utilised secondary definitions that are defined locally in the “produce_submod.py” file. Custom matplotlib styles (e.g. custom.mplstyle) used are defined in the “style” subfolder.

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
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import BoundaryNorm

from .produce_submod import *


def fig_1(data_dir,fig_dir,datafile='fig01.csv',figname='fig01',fmt='pdf',dpi='figure'):
    """
    
    Python script used to create Figure 1 in Argles, et al., in review.
    
    Caption
    -------
        
    Fig 1. Adapted from Shugart et al., 2020 [48], schematic showing how vegetation processes scale across space and time. Grey boxes classify the biological and community scale. Yellow text describes leaf-based processes. Blue text describes endogenous processes. Red text describes exogenous disturbances and processes. The black dotted lines illustrate the current spatial and temporal ranges of Earth System Models.

    References
    ----------
    
    Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

    [48] Shugart H, Foster A, Wang B, Druckenbrod D, Ma J, Lerdau M, et al. Gap models across micro-to mega-scales of time and space: examples of Tansley’s ecosystem concept. Forest Ecosystems. 2020;7(1):1–18.

    Parameters
    ----------
    data_dir : str
        Where the figure data is stored.
    fig_dir : str
        Where the figure is saved.
    datafile : str, optional
        The filename of the figure data. For the correct file formatting see below or in "data_dir/readme.txt" within the fig01 section. The default is 'fig01.csv'.
    figname : str, optional
        The name of the file saved in the figure directory. The default is 'fig01'.
    fmt : str, optional
        The type of formatting for the output figure file. Please refer to matplotlib.pyplot.savefig for the allowed file formats. The default is 'pdf'.
    dpi : float or 'figure', optional
        dots per inch, only applicable to relevant formats (e.g. fmt='png'). The default is 'figure'.

    """
    
    ### Data
    # Here we read in the figures data.
    #
    # datafile Format
    # ---------------
    #  For the correct formatting, datafile should be a csv file with the column headers: "feature name", "space center (m2)", "time center (yr)", "min space (m2)", "max space (m2)", "min time (yr)", "max time (yr)", "plot feature", "text pos", "text rotation", "zorder"
    #
    #  In a broad/rough sense; "space center (m2)", "time center (yr)", "min space (m2)", "max space (m2)", "min time (yr)", "max time (yr)" are meant to represent the spatial and temporal scales of biogeographical and ecological categorisation and processes.
    # 
    #  "plot feature" == "Rectangle" is used to describe biological and community spatial-temporal scales.
    #  "plot feature" == "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)" is used to describe a processes. 
    #  "plot feature" == "Line" this in the default figure refers to a timeframe (e.g. hour, day, year) or a geographic (e.g. Amazon Rainforest) extent.
    #  "plot feature" == "Line (model)" is used to describe a modelling (e.g. CMIP) extent.
    # 
    # Column Definitions
    # ------------------
    # feature name : str
    #	   The annotation attached to the object plotted.
    # space center (m2) : float
    #   The central space value of the annotation. Only used when column value at "plot feature" is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
    # time center (yr) : float
    #   The central time value of the annotation. Only used when column value at "plot feature" is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
    # min space (m2) : float
    #   The minimum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
    # max space (m2) : float
    #   The maximum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
    # min time (yr) : float
    #   The minimum value of time used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
    # max time (yr) : float
    #   The maximum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
    # plot feature : str
    #   The type of object placed on the plot. Can be equal to "Rectangle", "Text (Leaf)", "Text (endogenous)", "Text (exogenous)", "Line" or "Line (model)".
    # text_pos : str
    #   Describes where to place a string relative to rectangle or line on plot. First word is vertical placement ("top", "center", or "bottom"), second word is horizontal placement ("left", "center", or "right), e.g. "top left" (there must be a space). Not employed when “plot feature” is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
    # text_rotation : float
    #   Describes the angle of annotation relative to positive x-axis. (degrees) 
    # zorder : int
    #   The drawing order of plotted objects (e.g. zorder = 1, will place object above zorder = 0 but beneath zorder > 1).
    
    veg_scale_df = pd.read_csv(data_dir+'/'+datafile)
    num_objects = len(veg_scale_df)
 
    # Other values used in plot.
    temporal_lims = np.array([1.0/365.25*0.5/24.0,300e6]) # yaxis limit
    spatial_lims = np.array([1e-2,5e14]) # xaxis limit (corresponds to earth Surface Area)
    
    
    ### Figure
    
    # Pick style for figure (stored locally)
    plot_style('custom.mplstyle',global_style=False)
    
    gs = gridspec.GridSpec(1,1) # Define number of subplots
    fig = plt.figure()    
    ax = fig.add_subplot(gs[:,:])
 
    # Set axis limits, scale, labels, and titles
    ax.set_ylim(temporal_lims)
    ax.set_xlim(spatial_lims)   
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'Temporal Scale $\left(\mathrm{yr}\right)$',fontsize=10)
    ax.set_xlabel(r'Spatial Scale $\left(\mathrm{m^2}\right)$',fontsize=10)
    ax.set_title('Dynamic Vegetation and Scale',fontsize=10)
    
    # Define colors and fontsizes used for plotted objects
    rectangle_edgecolor = [0,0,0]
    rectangle_facecolor = [0,0,0,0.1]
    endogenous_color =  [0/255, 68/255, 136/255]
    exogenous_color = [187/255, 85/255, 102/255]
    leaf_color = [221/255, 170/255, 51/255]    
    rectangle_fontsize = 10   
    text_fontsize = 8
    
    # For loop runs through each of the saved objects in data file, here we plot each scale and process across space (x-axis) and time (y-axis).
    for n in range(0,num_objects):
    
        obj = veg_scale_df.iloc[n,:]            
        obj_name = obj['feature name']
        obj_zorder = obj['zorder']
        text_rotation = obj['text rotation']
        
        # For specific annotations we want to use two lines instead of one to save on space (for some reason "\n" is not properly used in ax.annotation() in via pandas).
        if obj_name == 'Stomatal Diurnal Processes':
            
            obj_name = 'Stomatal\nDiurnal Processes'
        
        elif obj_name == 'Gap-Phase Competition':
            
            obj_name = 'Gap-Phase\nCompetition'

        elif  obj_name == 'Glacial-Interglacial Climate Cycles':
            
            obj_name = 'Glacial-\nInterglacial\nClimate Cycles'
            
        #  "Rectangle" is used to describe biological and community spatial-temporal scales.
        if obj['plot feature'] == 'Rectangle':
        
            min_space = obj['min space (m2)']
            max_space = obj['max space (m2)']
            min_time = obj['min time (yr)']
            max_time = obj['max time (yr)']
            width = max_space - min_space
            height = max_time - min_time
            lower_right = (min_space,min_time)
            text_pos, text_ha, text_va = \
                interpret_text_pos(obj['text pos'],\
                                   (min_space,max_space),(min_time,max_time),\
                                   scale='log')
            rectangle = \
                patches.Rectangle(lower_right,width,height,linewidth=0.25,\
                                  edgecolor=rectangle_edgecolor,\
                                  facecolor=rectangle_facecolor,\
                                  zorder=obj_zorder)
            ax.add_patch(rectangle)
            ax.annotate(obj_name,xy=text_pos,ha=text_ha,va=text_va,\
                        fontsize=rectangle_fontsize,zorder=obj_zorder,\
                        color=rectangle_edgecolor,rotation=text_rotation,\
                        annotation_clip=False)  
    
        # "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)" is used to describe a processes. 
        elif obj['plot feature'] == 'Text (Leaf)':
            
            space_center = obj['space center (m2)']
            time_center = obj['time center (yr)']
            ax.annotate(obj_name,xy=(space_center,time_center),\
                        ha='center',va='center',color=leaf_color,\
                        fontsize=text_fontsize,zorder=obj_zorder,\
                        weight="bold",rotation=text_rotation,\
                        annotation_clip=False)
            
        elif obj['plot feature'] == 'Text (endogenous)':
            
            space_center = obj['space center (m2)']
            time_center = obj['time center (yr)']
            ax.annotate(obj_name,xy=(space_center,time_center),\
                        ha='center',va='center',color=endogenous_color,\
                        fontsize=text_fontsize,zorder=obj_zorder,\
                        weight="bold",rotation=text_rotation,\
                        annotation_clip=False)
        
        elif obj['plot feature'] == 'Text (exogenous)':
            
            space_center = obj['space center (m2)']
            time_center = obj['time center (yr)']
            ax.annotate(obj_name,xy=(space_center,time_center),\
                        ha='center',va='center',color=exogenous_color,\
                        fontsize=text_fontsize,zorder=obj_zorder,\
                        weight="bold",rotation=text_rotation,\
                        annotation_clip=False)

        # "Line"  refers to a timeframe (e.g. hour, day, year) or a geographic (e.g. Amazon Rainforest) extent.
        elif obj['plot feature'] == 'Line':
                        
            min_space = obj['min space (m2)']
            max_space = obj['max space (m2)']
            min_time = obj['min time (yr)']
            max_time = obj['max time (yr)']                    
            text_pos, text_ha, text_va = \
                interpret_text_pos(obj['text pos'],\
                                   (min_space,max_space),(min_time,max_time),\
                                   scale='log')
                            
            ax.plot([min_space,max_space],[min_time,max_time],\
                    color='black',lw=0.25)                            
            ax.annotate(obj_name,xy=text_pos,ha=text_ha,va=text_va,\
                        fontsize=text_fontsize,zorder=obj_zorder,\
                        color=rectangle_edgecolor,rotation=text_rotation,\
                        annotation_clip=False)
        
        # "Line (model)" is used to describe a modelling (e.g. CMIP) extent.
        elif obj['plot feature'] == 'Line (model)':
                        
            min_space = obj['min space (m2)']
            max_space = obj['max space (m2)']
            min_time = obj['min time (yr)']
            max_time = obj['max time (yr)']                    
            text_pos, text_ha, text_va = \
                interpret_text_pos(obj['text pos'],\
                                   (min_space,max_space),(min_time,max_time),\
                                   scale='log')
                            
            ax.plot([min_space,max_space],[min_time,max_time],\
                    color='black',ls=':')                            
            ax.annotate(obj_name,xy=text_pos,ha=text_ha,va=text_va,\
                        fontsize=text_fontsize,zorder=obj_zorder,\
                        color=rectangle_edgecolor,rotation=text_rotation,\
                        annotation_clip=False)    
                
    # Adjust subplot to maximise space in saved figure
    plt.subplots_adjust(left=0.075, bottom=0.105, right=0.9995, top=0.95)
    
    # Save figure
    plt.savefig(fig_dir+'/'+figname+'.'+fmt,dpi=dpi)
    
    return

def fig_2(data_dir,fig_dir,datafile='fig02.pickle',figname='fig02',fmt='pdf',dpi='figure'):
    """
    
    Python script used to create Figure 2 in Argles, et al., in review.
    
    Caption
    -------
        
    Fig 2. A rudimentary simulation of a forest from bare soil for 600 years, using the continuous patch eq (2) along with a toy gap model using a fixed 50x50 grid of patches on an arbitrary scale, and for a patch disturbance rate of  = 0:01yr􀀀1. Panel (a) shows how the age-distribution of patches changes through time for both the gap and statistical model. (b) shows the time-evolving spatial age map for the 50x50 grid.

    References
    ----------
    
    Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

    Shugart H, Foster A, Wang B, Druckenbrod D, Ma J, Lerdau M, et al. Gap models across micro-to mega-scales of time and space: examples of Tansley’s ecosystem concept. Forest Ecosystems. 2020;7(1):1–18.

    Parameters
    ----------
    data_dir : str
        Where the figure data is stored.
    fig_dir : str
        Where the figure is saved.
    datafile : str, optional
        The filename of the figure data. For the correct file formatting see below or in "data_dir/readme.txt" within the fig02 section. The default is 'fig02.pickle'.
    figname : str, optional
        The name of the file saved in the figure directory. The default is 'fig02'.
    fmt : str, optional
        The type of formatting for the output figure file. Please refer to matplotlib.pyplot.savefig for the allowed file formats. The default is 'pdf'.
    dpi : float or 'figure', optional
        dots per inch, only applicable to relevant formats (e.g. fmt='png'). The default is 'figure'.

    """

    ### Data
    # Here we read in the figures data.
    #
    # The data for this figure is created by running "generate_fig_2_patch_data" in the "methods_mod.py" file in the local "routines" module.    
    #
    # datafile Format
    # ---------------
    
    # For the correct formatting, datafile should be a pickle file with the keys: "t_slices", "k_sclices", "grid_res", "t", "dt", "loss_rate", "age_VF", "patch_age_pdf_VF", "num_patches_VF", "patch_age_stochastic_grid", "age_bins_S", "age_bin_centers_S", "age_bin_widths_S", "patch_age_pdf_S" 
    #
    # Key Definitions
    # ------------------
    # t_slices :  array_like[size=(N,),dtype=float], optional
    #    Where we want to save the spatial and age-structure time-slices for fig 2.  
    # k_sclices : array_like[size=(K_slice,),dtype=int],
    #    The k indicies of the time-slice for t. (-)
    # grid_res :  array_like[size=(2,),dtype=int], 
    #    The grid-resolution of the random run, the total number of patches simulated is the product of the elements.
    # t : array_like(size=(K,),dtype=float)
    #    The current time of the simulation. (yr)
    # dt : float
    #    The timestep of the simulation. (yr)
    # loss_rate : float
    #    The rate of patch disturbance. (/yr)
    # age_VF : array_like(size(K,),dtype=float)
    #    The Von Foerster patch age dimension. (yr)
    # patch_age_pdf_VF : array_like[size=(K_slice,K),dtype=float]
    #    The Von Foerster patch age probablity density function across the time slices. (/yr)
    # num_patches_VF : array_like(size(K,),dtype=int)
    #    The number of ages with positive pdf over simulation time. (-)
    # patch_age_stochastic_grid : array_like[size=(K_slice,I,J),dtype=float]
    #    The stochastic gridded patches across the time slices. (/yr) 
    # age_bins_S : array_like[size(K_slice,),dtype=object]
    #    The bin boundaries of the histogram across the time sclices.
    # age_bin_centers_S : array_like[size(K_slice,),dtype=object]
    #    The bin centers of the histogram across the time slices. (yr)
    # age_bin_widths_S : array_like[size(K_slice,),dtype=object] 
    #    The histogram bin width across the time slices. (yr)
    # patch_age_pdf_S : array_like[size(K_slice,),dtype=object]
    #    The stochastic grid patch age probablity density function of the histogram across the time slices. (/yr)

    with open(data_dir+'/'+datafile, 'rb') as handle:
        
        data_dict = pickle.load(handle)

    t_slices = data_dict['t_slices']
    k_sclices = data_dict['k_sclices']
    grid_res = data_dict['grid_res']
    I, J = grid_res[0], grid_res[1]
    t = data_dict['t'] 
    dt = data_dict['dt']
    loss_rate = data_dict['loss_rate']
    K = len(t)
    K_slice = len(t_slices)
    age_VF = data_dict['age_VF'] 
    patch_age_pdf_VF = data_dict['patch_age_pdf_VF']
    num_patches_VF = data_dict['num_patches_VF']
            
    patch_age_stochastic_grid = data_dict['patch_age_stochastic_grid']
    age_bins_S = data_dict['age_bins_S']
    age_bin_centers_S = data_dict['age_bin_centers_S']
    age_bin_widths_S = data_dict['age_bin_widths_S']
    patch_age_pdf_S = data_dict['patch_age_pdf_S']  
    
    
    roman_num = ['i','ii','iii','iv','v','vi']
    
    ### Figure
    
    # Pick style for figure (stored locally)
    plot_style('custom.mplstyle',global_style=False)
    gs = gridspec.GridSpec(4,K_slice,height_ratios=(1.0,0.55,1.0,0.1),\
                           hspace=0.1,wspace=0.1) # Define number of subplots
    fig = plt.figure()
    
    # Place figure title across all subplots
    ax = fig.add_subplot(gs[0,:])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('a) Patch Age Structure',pad=15)
    ax.set_xlabel(r'Patch Age $\left(\mathrm{yr}\right)$',labelpad=17.5)
    ax.set_ylabel(r'Patch Density $\left(\mathrm{yr^{-1}}\right)$',labelpad=25)
        
    for k in range(0,K_slice):
        
        ax = fig.add_subplot(gs[0,k])
        ax.set_ylim(1e-5,0.2)
        ax.set_xlim(0,500)
        ax.set_yscale('log')
        ax.set_title('a.'+str(roman_num[k])+') '+str(int(t_slices[k]))+r' yr')
        
        if k != 0:
            
            ax.set_yticks([])
    
        ax.plot(age_VF,patch_age_pdf_VF[k],lw=2.0,color='#004488',\
                label='Statistical Model')
        ax.bar(age_bin_centers_S[k],patch_age_pdf_S[k],age_bin_widths_S[k],\
               color='#BB5566',label='Stochastic Model')
            
        if k == 0:
            
            ax.annotate(r'$\Delta t=$ '+str(round(dt,2))+\
                        r' $\mathrm{yr}$'+'\n'+'  $\lambda=$ '+str(loss_rate)+\
                        r' $\mathrm{yr^{-1}}$',xy=(300,1e-3),\
                        ha='center',va='center')
            ax.legend(loc='upper right')
        
    
    cmap = plt.get_cmap('BrBG')
    norm = BoundaryNorm(np.arange(0,40,step=1), ncolors=cmap.N,extend='max')

    ax = fig.add_subplot(gs[2,:])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('b) Stochastic Model Canopy Output',pad=15)

    for k in range(0,K_slice):
        
        ax = fig.add_subplot(gs[2,k])
        ax.set_title('b.'+str(roman_num[k])+') '+str(int(t_slices[k]))+r' yr')
        ax.set_yticks([])
        ax.set_xticks([])
    
        c = ax.pcolormesh(patch_age_stochastic_grid[k,:,:], cmap=cmap, norm=norm)
        
    
    ax = fig.add_subplot(gs[3,:])
    plt.colorbar(c,cax=ax,orientation='horizontal')
    ax.set_xlabel('Patch Age')
    
    # Adjust subplot to maximise space in saved figure
    plt.subplots_adjust(left=0.075, bottom=0.085, right=0.995, top=0.925)
    
    # Save figure
    plt.savefig(fig_dir+'/'+figname+'.'+fmt,dpi=dpi)
    
    
    return

def fig_4(data_dir,fig_dir,datafile='fig04.csv',figname='fig04',fmt='pdf',dpi='figure'):
    """
    
    Python script used to create Figure 4 in Argles, et al., in review.
    
    Caption
    -------
        
    Fig 4. Dynamic Global Vegetation Model (DGVM) trends through time. Here models are categorized into three types: Average, Intermediate, and Complex. (1) Simple - plant size or age demographics are not included beyond the "average" case. (2) Intermediate - models include the representation of size and/or age within cohorts. (3) Complex - individual based models.

    References
    ----------
    
    Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

    Parameters
    ----------
    data_dir : str
        Where the figure data is stored.
    fig_dir : str
        Where the figure is saved.
    datafile : str, optional
        The filename of the figure data. For the correct file formatting see below or in "data_dir/readme.txt" within the fig04 section. The default is 'fig04.csv'.
    figname : str, optional
        The name of the file saved in the figure directory. The default is 'fig04'.
    fmt : str, optional
        The type of formatting for the output figure file. Please refer to matplotlib.pyplot.savefig for the allowed file formats. The default is 'pdf'.
    dpi : float or 'figure', optional
        dots per inch, only applicable to relevant formats (e.g. fmt='png'). The default is 'figure'.

    """
    
    ### Data
    # Here we read in the figures data.
    #
    # datafile Format
    # ---------------
    #  For the correct formatting, datafile should be a csv file with the column headers: "Model", "Method", "Publication Year", "Citations"
    # 
    # Column Definitions
    # ------------------
    # Model : str
    #	   The DGVM model name or acronym.
    # Method : str
    #   The generalised catergorisation of the DGVM.
    # Publication Year : float
    #   The year of model paper publication.
    # Citations : str
    #   The citations that are associated with the paper, for full citations see references in Argles, et al., in review.
        
    dgvm_data = pd.read_csv(data_dir+'/'+datafile)
    year = dgvm_data['Publication Year']
    num_dgvms = len(year)
    
    ### Figure
    
    # Pick style for figure (stored locally)
    plot_style('custom.mplstyle',global_style=False)
    gs = gridspec.GridSpec(2,2,height_ratios=(0.4,1),hspace=0.02) # Define number of subplots
    fig = plt.figure()
    ax = fig.add_subplot(gs[:,:])    
    
    # Place figure title across all subplots
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('DGVMs over time',pad=15)
        
    # Set axis limits, scale, labels, and titles
    
    ax = fig.add_subplot(gs[0,:])
    xlim_min = 1995 # year before first publication
    xlim_max = np.nanmax(year)+0.5
    ylim_min = 0
    ylim_max = len(year[year==year])+0.5+2
    ax.set_xlim(xlim_min,xlim_max)
    ax.set_ylim(0.8,4.2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    years = np.linspace(xlim_min,2020,2020-xlim_min+1)
    
    # Assign counts to each method
    num_years = len(years)    
    total_count = np.zeros(num_years,dtype=int)
    individual_count = np.zeros(num_years,dtype=int)
    average_count = np.zeros(num_years,dtype=int)
    cohort_count = np.zeros(num_years,dtype=int)
    cohort_1d_count = np.zeros(num_years,dtype=int)
    
    # Define colors for different DGVM methods
    current_year_count = 1
    last_year = np.nan
    colors = ['#004488','#BB5566','#DDAA33']
    
    # For loop runs each DGVM saved in the datafile counting for each method defined. We also plot on top panel the name and colorcoded point for each DGVM.
    
    for dgvm in range(0,num_dgvms):
            
        method = dgvm_data['Method'].iloc[dgvm]        
    
        if year[dgvm] == year[dgvm]:
        
           # Default marker panel values, DGVMs with no defined method will be proscribed a black marker on top panel
           avg_ind_mixed = False
           avg_coh_mixed = False
           marker = 'o'
           marker_color = 'black'
           marker_fillstyle = 'full'
           
           # If defined method (not empty value) we count the total for that year and for the specific DGVM method 
           if method == method:
               
               total_count[years==year[dgvm]] = total_count[years==year[dgvm]] + 1
                 
               if 'Averaged' in method:
                       
                   average_count[years==year[dgvm]] = average_count[years==year[dgvm]] + 1
                   marker = 'o'
                   marker_color = colors[1]
                   
                   if 'Individual' in method:
                       
                       total_count[years==year[dgvm]] = total_count[years==year[dgvm]] + 1
                       individual_count[years==year[dgvm]] = individual_count[years==year[dgvm]] + 1
                       avg_ind_mixed = True 
               
                   if 'Cohort' in method:
                       
                       total_count[years==year[dgvm]] = total_count[years==year[dgvm]] + 1
                       cohort_count[years==year[dgvm]] = cohort_count[years==year[dgvm]] + 1
                       avg_coh_mixed = True
                       
                
               elif 'Individual' in method:
                   
                   marker = 'o'
                   marker_color = colors[0]
                   individual_count[years==year[dgvm]] = individual_count[years==year[dgvm]] + 1
               
               elif 'Cohort' in method:
                   
                   cohort_count[years==year[dgvm]] = cohort_count[years==year[dgvm]] + 1 
                   marker = 'o'
                   marker_color = colors[2]
                   
                   if 'and' not in method:
                       marker_fillstyle = 'left'  
                       cohort_1d_count[years==year[dgvm]] = cohort_1d_count[years==year[dgvm]] + 1
               
               elif 'Climate' in method:
                    
                    marker = '*'
                    
           if year[dgvm] == last_year:
                
               current_year_count = current_year_count + 1
                
           else:
                
               current_year_count = 1
                   
           #  Plot DGVM scatter point and label. Height on y-axis is dependent on the number of DGVMs already there (current_year_count).
           if (avg_ind_mixed == False and avg_coh_mixed == False):
    
               ax.plot(year[dgvm],current_year_count,zorder=1,marker=marker,markeredgecolor=marker_color,markerfacecolor=marker_color,fillstyle=marker_fillstyle)
               ax.annotate(dgvm_data['Model'].iloc[dgvm],xy=(year[dgvm]+0.2,current_year_count+0.04),rotation=45,zorder=2,fontsize=8)
               last_year = year[dgvm]
               
            # Some DGVMs employ multiple modes (such as LPG-GUESS) here we want to place in multiple method caterogies, we also modify the market to be different colors
           elif (avg_ind_mixed == True and avg_coh_mixed == True):
               
               ax.plot(year[dgvm],current_year_count,zorder=3,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[1],fillstyle='left')
               ax.plot(year[dgvm],current_year_count,zorder=2,markersize=7,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[1],fillstyle='left')
               ax.plot(year[dgvm],current_year_count,zorder=1,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[0],fillstyle='top')
               ax.plot(year[dgvm],current_year_count,zorder=0,markersize=7,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[0],fillstyle='top')
               ax.plot(year[dgvm],current_year_count,zorder=1,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[2],fillstyle='bottom')
               ax.plot(year[dgvm],current_year_count,zorder=0,markersize=7,marker=marker,markeredgecolor=[0,0,0,0],markerfacecolor=colors[2],fillstyle='bottom')
               ax.annotate(dgvm_data['Model'].iloc[dgvm],xy=(year[dgvm]+0.2,current_year_count+0.04),rotation=45,zorder=2,fontsize=8)
               last_year = year[dgvm]
               
    # Find the cummulative total for each method across year.
    total_count = np.cumsum(total_count)
    individual_count = np.cumsum(individual_count)
    average_count = np.cumsum(average_count)
    cohort_count = np.cumsum(cohort_count)
    cohort_1d_count = np.cumsum(cohort_1d_count)      
    
    # Bottom panel is for the general trend of each method catergory across year. 
    ax = fig.add_subplot(gs[1,:])
    ax.set_xlim(xlim_min,xlim_max)
    ax.set_ylim(ylim_min,20)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.plot(years,individual_count,zorder=1,linewidth=1.5,label='Individual Based')
    ax.plot(years,average_count,zorder=2,linewidth=1.5,label='Average Based')
    ax.plot(years,cohort_count,zorder=3,linewidth=1.5,label='Cohort Based')
    ax.plot(years,cohort_1d_count,color=colors[2],zorder=3,linewidth=1.5,ls='--',label='Cohort 1D Subtotal')
    ax.legend(loc='center left')
    ax.spines['top'].set_visible(True)
    ax.tick_params(top='on',which='both')
    
    ax.set_xlabel('Publication Year')
    ax.set_ylabel('Count')
        
    
    # Additionally we want to plot a simple annotation detailing the "complexity" of each method.
    ax.annotate('', xy=(1998,15.5), xytext=(2006,15.5), arrowprops=dict(arrowstyle='<->',linewidth=1.5))
    diff = 2006 - 1998.35
    ax.annotate('"Complex"',xy=(1996.75,14.5))
    ax.annotate('"Simple"',xy=(2005.35,14.5))
    ax.annotate('"Intermediate"',xy=(2000.5,14.5))
    ax.plot([1998.35],[17.25],marker='o',markeredgecolor=colors[0],markerfacecolor=colors[0])
    ax.annotate('Individual',xy=(1997.45,17.75))
    ax.plot([1998.35+diff*1.0/3.0],[17.25],marker='o',markeredgecolor=colors[2],markerfacecolor=colors[2])
    ax.annotate('Cohort 2D',xy=(1998.35+diff*1.0/3.0-1,16))
    ax.plot([1998.35+diff*2.0/3.0],[17.25],marker='o',markeredgecolor=colors[2],markerfacecolor=colors[2],fillstyle='left')
    ax.annotate('Cohort 1D',xy=(1998.35+diff*2.0/3.0-1,17.75))
    ax.plot([2005.6],[17.25],marker='o',markeredgecolor=colors[1],markerfacecolor=colors[1])
    ax.annotate('Average',xy=(2005.65-0.8,16.))
    ax.set_yticks([0,5,10,15,20])
      
    # Adjust subplot to maximise space in saved figure
    plt.subplots_adjust(left=0.0545, bottom=0.0835, right=0.9625, top=0.925)
    
    # Save figure
    plt.savefig(fig_dir+'/'+figname+'.'+fmt,dpi=dpi)
    
    return


