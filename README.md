# Turtle-Labyrinth-Ecomorphology-and-Evolution-Data
Supplementary Data Files to the paper "Independent origin of large labyrinth size in turtles"

TO PERFORM ANALYSES ON TURTLE LABYRINTH DATA

1. Store Datasets S1–6: these contain data to be read by the R-scripts
	- Dataset S1: specimen names, ID, measurements, ecological classifications, other info for turtle data
	- Dataset S2: Landmark data for 168 turtle specimens as csv files (3D coordinate data)
	- Dataset S3: information about which landmarks are used as sliding semilandmarks
	- Dataset S4: information about landmark colours used in one R script for deformation plots
	- Dataset S5: phylogenetic tree based on cal3-calibration, used for main analyses
	- Dataset S6: phylogenetic tree based on mil-calibration, used for sensitivity analyses
2. Open Datasets S7–9 & 12 in R to run different types of analyses. These scripts contain commands to read Datasets S1–6.
	- Dataset S7: performs GPA & PCA, produces plots for Figure 1 of main text
	- Dataset S8: performs GPA & labyrinth shape regression analyses, produces Table 1 of main text
	- Dataset S9: performs GPA & labyrinth size regression analyses, produces Table 2 of main text
		- Dataset S10: output of S9, represents full version of table printed as Table 2 of main text
		- Dataset S11: as S10, but using different phylogenetic tree as a sensitivity analysis
	- Dataset S12: prints regression and residual plot, and performs ancestral state reconstructions; produces Figure 2 of main text
3. Note that Datasets S10–11 are an output list of the script Dataset S9


#####
TO PERFORM ANALYSES ON SMNIOTE LABYRINTH DATA

1. Store Datasets S13–15: these contain data to be read by the R-scripts
	- Dataset S13: specimen names, ID, measurements, ecological classifications, other info for amniote data
	- Dataset S14: Landmark data for 200 amniote specimens as csv files (3D coordinate data)
	- Dataset S15: information about which landmarks are used as sliding semilandmarks
2. Open Dataset S16 in R to run the amniote analysis. This script contains commands to read Datasets S13–15.
