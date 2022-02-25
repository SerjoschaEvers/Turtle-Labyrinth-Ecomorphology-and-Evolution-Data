# Turtle-Labyrinth-Ecomorphology-and-Evolution-Data
Supplementary Data Files to the paper "Independent origin of large labyrinth size in turtles"

#####

v1.1 Update with regard to v1.0: updated taxonomy for mammal species in files to match MorphoSource records of CT data

#####

TO PERFORM ANALYSES ON TURTLE LABYRINTH DATA

1. Store Datasets 1–6: these contain data to be read by the R-scripts
	- Dataset 1: specimen names, ID, measurements, ecological classifications, other info for turtle data
	- Dataset 2: Landmark data for 168 turtle specimens as csv files (3D coordinate data)
	- Dataset 3: information about which landmarks are used as sliding semilandmarks
	- Dataset 4: information about landmark colours used in one R script for deformation plots
	- Dataset 5: phylogenetic tree based on cal3-calibration, used for main analyses
	- Dataset 6: phylogenetic tree based on mil-calibration, used for sensitivity analyses
2. Open Datasets 7–9 & 12 in R to run different types of analyses. These scripts contain commands to read Datasets S1–6.
	- Dataset 7: performs GPA & PCA, produces plots for Figure 1 of main text
	- Dataset 8: performs GPA & labyrinth shape regression analyses, produces Table 1 of main text
	- Dataset 9: performs GPA & labyrinth size regression analyses, produces Table 2 of main text
		- Dataset 10: output of 9, represents full version of table printed as Table 2 of main text
		- Dataset 11: as S10, but using different phylogenetic tree as a sensitivity analysis
	- Dataset 12: prints regression and residual plot, and performs ancestral state reconstructions; produces Figure 2 of main text
3. Note that Datasets 10–11 are an output list of the script Dataset 9


#####
TO PERFORM ANALYSES ON SMNIOTE LABYRINTH DATA

1. Store Datasets 13–15: these contain data to be read by the R-scripts
	- Dataset 13: specimen names, ID, measurements, ecological classifications, other info for amniote data
	- Dataset 14: Landmark data for 200 amniote specimens as csv files (3D coordinate data)
	- Dataset 15: information about which landmarks are used as sliding semilandmarks
2. Open Dataset 16 in R to run the amniote analysis. This script contains commands to read Datasets 13–15.
