# code_Argles_PLOS_Climate_22
This contains the data and code to produce figures 1, 2 and 4 in Dynamic Global Vegetation Models: searching for the balance between demographic process representation and computational tractability (Argles et al., in review PLOS Climate.). 

# How to run
The "main.py" generates some and all figure data for 1 and 2 respectfully. The methods used to generate this data are defined in “routines/methods_mod.py” and “routines/methods_submod.py”. In the "data" subfolder we have put all the data required to reproduce the figures present in the review paper, along with a brief description of the data sources and assumptions in "data/readme.txt".  The plotting routines are defined in “routines/produce_mod.py” and “routines/produce_submod.py”, with the resultant figures saved in “figures/”. Custom matplotlib files are defined in “routines/styles/”.

# Figure Captions

Fig 1. Adapted from Shugart et al., 2020 [48], schematic showing how vegetation processes scale across space and time. Grey boxes classify the biological and community scale. Yellow text describes leaf-based processes. Blue text describes endogenous processes. Red text describes exogenous disturbances and processes. The black dotted lines illustrate the current spatial and temporal ranges of Earth System Models.

Fig 2. A rudimentary simulation of a forest from bare soil for 600 years, using the continuous patch eq (2) along with a toy gap model using a fixed 50x50 grid of patches on an arbitrary scale, and for a patch disturbance rate of loss_rate = 0:01/yr. Panel (a) shows how the age-distribution of patches changes through time for both the gap and statistical model. (b) shows the time-evolving spatial age map for the 50x50 grid.

Fig 4. Dynamic Global Vegetation Model (DGVM) trends through time. Here models are categorized into three types: Average, Intermediate, and Complex. (1) Simple - plant size or age demographics are not included beyond the "average" case. (2) Intermediate - models include the representation of size and/or age within cohorts. (3) Complex - individual based models.

# References
   
   
Argles, A.P.K., Moore, J.R., and Cox, P.M. Dynamic Global Vegetation Models: searching or the balance between demographic process representation and computational tractability. In review (PLOS Climate).

[48] Shugart H, Foster A, Wang B, Druckenbrod D, Ma J, Lerdau M, et al. Gap models across micro-to mega-scales of time and space: examples of Tansley’s ecosystem concept. Forest Ecosystems. 2020;7(1):1–18.
