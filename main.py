# -*- coding: utf-8 -*-
"""
Created on 22/04/2022

This Python script is used to generate figures seen in the literature review: Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability (Argles, et al., in review PLOS Climate).

The script saves the following figures in the 'figures' subfolder:
    
    Fig 1. Adapted from Shugart et al., 2020 [48], schematic showing how vegetation processes scale across space and time. Grey boxes classify the biological and community scale. Yellow text describes leaf-based processes. Blue text describes endogenous processes. Red text describes exogenous disturbances and processes. The black dotted lines illustrate the current spatial and temporal ranges of Earth System Models.

    Fig 2. A rudimentary simulation of a forest from bare soil for 600 years, using the continuous patch eq (2) along with a toy gap model using a fixed 50x50 grid of patches on an arbitrary scale, and for a patch disturbance rate of  = 0:01yr􀀀1. Panel (a) shows how the age-distribution of patches changes through time for both the gap and statistical model. (b) shows the time-evolving spatial age map for the 50x50 grid.

    Fig 4. Dynamic Global Vegetation Model (DGVM) trends through time. Here models are categorized into three types: Average, Intermediate, and Complex. (1) Simple - plant size or age demographics are not included beyond the "average" case. (2) Intermediate - models include the representation of size and/or age within cohorts. (3) Complex - individual based models.

Fig 3. from the paper is saved independently of the script as the results have been saved in other publications (Moore et al., 2018, 2020). For accesses to the dataset used in Figure 3 see the USDA (Oswalt, et al., 2010) and RAINFOR (Peacock, et al., 2007) datasets, and contact Dr Jonathan Moore for access to the code (j.moore3@exeter.ac.uk).

    Fig 3. Adapted from the work of Moore et al., 2018, 2020 [142, 143]. Best fits of Demographic Equilibrium Theory (DET), eq (5), for mu_0 (assuming z = d and phi = 1/3) against observations of number density versus diameter. Panel (a) shows the DET fit to North American trees [142] observations of trunk diameter (diameter at breast height) from the USDA Forest Service FIA program [144]. Panel (b) presents DET fits to RAINFOR measurements of basal diameter (b) [145]. mu_1 is the normalised metric of the turnover parameter (mu_0 at d0 = 1 cm) describing the shape of the Weibull distribution.

Additionally data used to make figures 1,2 and 4 are saved in the 'data' subfolder.

For queries please contact Dr Arthur Argles (A.Argles2@exeter.ac.uk) 

References
----------

Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

[48] Shugart H, Foster A, Wang B, Druckenbrod D, Ma J, Lerdau M, et al. Gap models across micro-to mega-scales of time and space: examples of Tansley’s ecosystem concept. Forest Ecosystems. 2020;7(1):1–18.

[142] Moore JR, Zhu K, Huntingford C, Cox PM. Equilibrium forest demography explains the distribution of tree sizes across North America. Environmental Research Letters. 2018;13(8):084019.

[143] Moore JR, Argles APK, Zhu K, Huntingford C, Cox PM. Validation of demographic equilibrium theory against tree-size distributions and biomass density in Amazonia. Biogeosciences. 2020;17(4):1013{1032. doi:10.5194/bg-17-1013-2020.
                                                                                                        
[144] Oswalt SN, Smith WB, Miles PD, Pugh SA. Forest Resources of the United States, 2012: a technical document supporting the Forest Service 2010 update of the RPA Assessment. Gen Tech Rep Wo-91 Washington, Dc: Us Department of Agriculture, Forest Service, Washington Office 218 P. 2014;91.

[145] Peacock J, Baker T, Lewis S, Lopez-Gonzalez G, Phillips O. The RAINFOR database: monitoring forest biomass and dynamics. Journal of Vegetation Science. 2007;18(4):535{542.

"""

# Credits
__author__ = "Arthur Argles, Jonathan Moore, and Peter Cox"
__credits__ = ["Arthur Argles","Jonathan Moore","Peter Cox"]
__license__ = "CCBY4.0"
__maintainer__ = "Arthur Argles"
__email__ = "A.Argles2@exeter.ac.uk"


import os
import numpy as np

# Local Import
import routines


if __name__ == '__main__':
    
    ### Directories
    # For the code directory structure we have the subfolders (other than __pycache__): 
    # 
    # 'code_Argles_PLOS_Climate_22' (1) 
    #    |
    #    |------------------------|------------------------|
    # 'routines' (2)           'data' (3)             'figures' (4)
    #    |
    # 'styles' (5)
    #
    #  (1) Where we have the main python script (main.py) to generate paper results
    #  (2) The routines used in the main script to simulate data (methods_mod.py) and produce figures (produce_mod.py)
    #  (3) Where the data utilised in producing the figures is stored.
    #  (4) Where the figures are saved.
    #  (5) Where any custom matplotlib styles (e.g. custom.mplstyle) are stored.  
    main_dir = os.getcwd()
    data_dir = os.path.join(main_dir,'data')
    fig_dir = os.path.join(main_dir,'figures')

    ### Data
    # Figure Data is stored in data_dir
    # Here we use demonstrate the methods used in generating certain figure results.
    #  
    # Other sources of data are provided in 'data/datainfo.txt', not all data is
    # generated here in methods_mod.py.
    #
    # Uncomment For function description
    help(routines.generate_fig_1_gridbox_area)
    help(routines.generate_fig_2_patch_data)    
    #
    # lon_res, lat_res is assumed to follow the N96 grid of the UK MO UM
    # lon_res = 1.88
    # lat_res = 1.25
    # routines.generate_fig_1_gridbox_area(data_dir,datafile='fig01.csv',\
    #                                      lon_res=lon_res,lat_res=lat_res)
    # Toy Patch Simulation Parameter Values        
    #
    # These values roughly follow the same time-step as the ED DGVM model (Moorcroft et al., 2001)
    t_end = 600.
    dt = 1.0/12.0
    loss_rate = 0.01
    t_slices = np.array([5.,50.,250.,500.])
    grid_res = [50,50]
    routines.generate_fig_2_patch_data(t_end,dt,loss_rate,data_dir,\
                                        datafile='fig02.pickle',\
                                        include_VF=True,include_grid=True,\
                                        grid_res=grid_res,t_slices=t_slices)
        
    ### Figures
    # Here we save figures 1, 2 and 4 in (Argles, et al. in review). 
    # Figure 3 was made by Dr Jonathan Moore (j.moore3@exeter.ac.uk) from results Moore et al., 2018,20 and is saved independently from this script. 
    #
    help(routines.fig_1)  
    help(routines.fig_2)  
    help(routines.fig_4)
    routines.fig_1(data_dir,fig_dir,fmt='pdf')
    routines.fig_2(data_dir,fig_dir,fmt='pdf')
    routines.fig_4(data_dir,fig_dir,fmt='pdf')