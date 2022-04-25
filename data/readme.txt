--------
Figure 1
--------

datafile Format
---------------
For the correct formatting, datafile should be a csv file with the column headers: "feature name", "space center (m2)", "time center (yr)", "min space (m2)", "max space (m2)", "min time (yr)", "max time (yr)", "plot feature", "text pos", "text rotation", "zorder"
 
In a broad/rough sense; "space center (m2)", "time center (yr)", "min space (m2)", "max space (m2)", "min time (yr)", "max time (yr)" are meant to represent the spatial and temporal scales of biogeographical and ecological categorisation and processes.
     
"plot feature" == "Rectangle" is used to describe biological and community spatial-temporal scales.
"plot feature" == "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)" is used to describe a processes. 
"plot feature" == "Line" this in the default figure refers to a timeframe (e.g. hour, day, year) or a geographic (e.g. Amazon Rainforest) extent.
"plot feature" == "Line (model)" is used to describe a modelling (e.g. CMIP) extent.
     
Column Definitions
------------------
feature name : str
	The annotation attached to the object plotted.
space center (m2) : float
      The central space value of the annotation. Only used when column value at "plot feature" is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
time center (yr) : float
	The central time value of the annotation. Only used when column value at "plot feature" is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
min space (m2) : float
	The minimum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
max space (m2) : float
	The maximum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
min time (yr) : float
      The minimum value of time used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
max time (yr) : float
	The maximum value of space used when plotting a Line or Rectangle. Only used when column value at "plot feature" is equal to "Rectangle", "Line", or "Line (model)" .
plot feature : str
	The type of object placed on the plot. Can be equal to "Rectangle", "Text (Leaf)", "Text (endogenous)", "Text (exogenous)", "Line" or "Line (model)".
text_pos : str
      Describes where to place a string relative to rectangle or line on plot. First word is vertical placement ("top", "center", or "bottom"), second word is horizontal placement ("left", "center", or "right), e.g. "top left" (there must be a space). Not employed when “plot feature” is equal to "Text (Leaf)", "Text (endogenous)" or "Text (exogenous)".
text_rotation : float
      Describes the angle of annotation relative to positive x-axis. (degrees) 
zorder : int
      The drawing order of plotted objects (e.g. zorder = 1, will place object above zorder = 0 but beneath zorder > 1).

Data Values
-------------

This figure is very much a "rough" outline of processes in ball-park figures. The main interpretation of the figure should be the "order" of processes in space and time. The data was interepreted from Figure 1 in Shugart H, et al., 2020 with a couple of more adaptations such as including leaf and ESM/DGVM target scales.

feature_name	
	
Leaf		min space (m2)	0.01		Note: A rough cut-off on the minimum surface area of a leaf, there are smaller leaves on certain plants.

Leaf		max space (m2)	2.82		Note: A rough cut-off on the maximum surface area of a leaf, again there are larger leaves on certain plants.

Leaf		min time  (yr)	5.70E-05	Note: half an hour in units yr. Half an hour is approx smallest timestep in Dynamic Global Vegetation Models.

Leaf		max time  (yr)	4		Note: 1/gamma_leaf_0 for tree PFTs in JULES.

Plant		min space (m2)	0.01		Note: Assume same min_space as Leaf, as some plants just have a few leaves.

Plant		max space (m2)	100		Note: A rough cut-off on the maximum crown area of a tree, mainly to save on space with the figure. Crown areas can exceed this (Martínez Cano, et al. 2019).

Stand, Landscape, Biome, Phenology,
Stomatal Diurnal Processes, Secondary Succession,
Primary Succession, Gap-Phase Competition, Recruitment
Evolution, Disturbance, Pathogens, Fire Regime, 
Human Activities, Climatic Fluctuations, Soil Development
Glacial-Interglacial Climate Cycles					Note: These processes are roughly taken from combining Figures 1.a-c in Shugart H, et al., 2020 there is a spread of processes in both space of time captured by the double headed arrows within this paper.

CMIP (250 yr) min/max time  (yr)	250	Note: CMIP historical and projections normally run between 1850 and 2100, i.e. 250 years.

Amazon Rainforest	min/max space (m2)	5.5E+12		Note: value from wikpedia (https://en.wikipedia.org/wiki/Amazon_rainforest accessed 22/04/2022) also found in various publications (e.g. Lorenz, et al., 2021)

Global Land Area	min/max space (m2)	1.49E+14		Note:	value from Pidwirny, M. (2006).

Gridbox Area Min	min/max space (m2)	1188927892		Note: Estimated closest grid-box area at a latitude 90°N assuming a spherical earth estimated in produce_mod.py (def generate_fig_1_gridbox_area). We assume N96 resolution (longitude resolution 1.88°, latitude resolution: 1.25°) from the UK MO UM model (Walters, et al., 2019).

Gridbox Area Max	min/max space (m2)	29117843884		Note: Estimated closest grid-box area at a latitude 0°N assuming a spherical earth estimated in produce_mod.py (def generate_fig_1_gridbox_area). We assume N96 resolution (longitude resolution 1.88°, latitude resolution: 1.25°) from the UK MO UM model (Walters, et al., 2019).


References
----------
Shugart, et al. Gap models across micro-to mega-scales of time and space: examples of Tansley’s ecosystem concept. Forest Ecosystems. 2020;7(1):1–18.

Martínez Cano, et al. Tropical tree height and crown allometries for the Barro Colorado Nature Monument, Panama: a comparison of alternative hierarchical models incorporating interspecific variation in relation to life history traits. Biogeosciences 16.4 (2019): 847-862.

Lorenz, et al. Deforestation hotspots, climate crisis, and the perfect scenario for the next epidemic: The Amazon time bomb. Science of the Total Environment 783 (2021): 147090.  

Pidwirny, M. Surface area of our planet covered by oceans and continents. University of British Columbia, Okanagan. 2006.

Walters, et al. The Met Office unified model global atmosphere 6.0/6.1 and JULES global land 6.0/6.1 configurations. Geoscientific Model Development. 2017;10(4):1487-1520.


--------
Figure 2
--------
 The data for this figure is created by running "generate_fig_2_patch_data" in the "methods_mod.py" file in the local "routines" module.    

 datafile Format
 ---------------
 For the correct formatting, datafile should be a pickle file with the keys: "t_slices", "k_sclices", "grid_res", "t", "dt", "loss_rate", "age_VF", "patch_age_pdf_VF", "num_patches_VF", "patch_age_stochastic_grid", "age_bins_S", "age_bin_centers_S", "age_bin_widths_S", "patch_age_pdf_S" 

 Key Definitions
 ------------------
 t_slices :  array_like[size=(N,),dtype=float], optional
	Where we want to save the spatial and age-structure time-slices for fig 2.  
 k_sclices : array_like[size=(K_slice,),dtype=int],
	The k indicies of the time-slice for t. (-)
 grid_res :  array_like[size=(2,),dtype=int], 
	The grid-resolution of the random run, the total number of patches simulated is the product of the elements.
 t : array_like(size=(K,),dtype=float)
	The current time of the simulation. (yr)
 dt : float
	The timestep of the simulation. (yr)
 loss_rate : float
	The rate of patch disturbance. (/yr)
 age_VF : array_like(size(K,),dtype=float)
	The Von Foerster patch age dimension. (yr)
 patch_age_pdf_VF : array_like[size=(K_slice,K),dtype=float]
	The Von Foerster patch age probablity density function across the time slices. (/yr)
 num_patches_VF : array_like(size(K,),dtype=int)
	The number of ages with positive pdf over simulation time. (-)
 patch_age_stochastic_grid : array_like[size=(K_slice,I,J),dtype=float]
	The stochastic gridded patches across the time slices. (/yr) 
 age_bins_S : array_like[size(K_slice,),dtype=object]
	The bin boundaries of the histogram across the time sclices.
 age_bin_centers_S : array_like[size(K_slice,),dtype=object]
	The bin centers of the histogram across the time slices. (yr)
 age_bin_widths_S : array_like[size(K_slice,),dtype=object] 
	The histogram bin width across the time slices. (yr)
 patch_age_pdf_S : array_like[size(K_slice,),dtype=object]
	The stochastic grid patch age probablity density function of the histogram across the time slices. (/yr)

--------
Figure 4
--------

datafile Format
---------------
For the correct formatting, datafile should be a csv file with the column headers: "Model", "Method", "Publication Year", "Citations"
 
Column Definitions
------------------
Model : str
	The DGVM model name or acronym.
Method : str
	The generalised catergorisation of the DGVM.
Publication Year : float
	The year of model paper publication.
Citations : str
	The citations that are associated with the paper.


References
----------

For full references + citations please refer to both Argles et al., in review and csv file.