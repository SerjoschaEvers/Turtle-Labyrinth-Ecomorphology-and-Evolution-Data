# Turtle-Labyrinth-Ecomorphology-and-Evolution-Data
Supplementary Data Files to the paper "Independent origin of large labyrinth size in turtles"

Read me

#####

v.3 Update with regard to v2.0: new script (Dataset 9) based on reviewer comments.

#####
TO PERFORM ANALYSES ON TURTLE LABYRINTH DATA

1. Store Datasets S1–6: these contain data to be read by the R-scripts
	- Dataset S1: specimen names, ID, measurements, ecological classifications, other info for turtle data
	- Dataset S2: Landmark data for 168 turtle specimens as csv files (3D coordinate data)
	- Dataset S3: information about which landmarks are used as sliding semilandmarks
	- Dataset S4: information about landmark colours used in one R script for deformation plots
	- Dataset S5: phylogenetic tree based on cal3-calibration, used for main analyses
	- Dataset S6: phylogenetic tree based on mil-calibration, used for sensitivity analyses
2. Open Datasets S7–10 & 13 in R to run different types of analyses. These scripts contain commands to read Datasets S1–6.
	- Dataset S7: performs GPA & PCA, produces plots for Figure 1 of main text
	- Dataset S8: performs GPA & labyrinth shape regression analyses, produces Table 1 of main text
	- Dataset S9: performs GPA & labyrinth shape regression analyses, produces Table 1 of main text
	- Dataset S10: performs GPA & size-corrected PCA, size-corrected shape regression analyses, produces Figure 2 of main text
		- Dataset S11: output of S10, represents full version of table printed as Table 2 of main text
		- Dataset S12: as S11, but using different phylogenetic tree as a sensitivity analysis
	- Dataset S12: prints regression and residual plot, and performs ancestral state reconstructions; produces Figure 3 of main text
3. Note that Datasets S11–12 are an output list of the script Dataset S10


#####
TO PERFORM ANALYSES ON AMNIOTE LABYRINTH DATA

1. Store Datasets S14–16: these contain data to be read by the R-scripts
	- Dataset S14: specimen names, ID, measurements, ecological classifications, other info for amniote data
	- Dataset S15: Landmark data for 200 amniote specimens as csv files (3D coordinate data)
	- Dataset S16: information about which landmarks are used as sliding semilandmarks
2. Open Dataset S17 in R to run the amniote analysis. This script contains commands to read Datasets S14–16. Produces Figure 4 of main text.
