# Turtle-Labyrinth-Ecomorphology-and-Evolution-Data
Supplementary Data Files to the paper "Independent origin of large labyrinth size in turtles"

Read me

#####

v.4 Update with regard to v3.0: new scripts (now Datasets 1, 6, 9) based on editorial comments that some tables were too long to be supplementary tables (Datasets 1, 6) and that we should move some of the supplementary methods into the actual methods (Dataset 9).

#####
Dataset 1: Contains lists of specimens with MorphoSource links to 3D models. CT scans are internally linked in MorphoSource. The spreadsheet also details the one vs. restricted download as currently implemented in MorphoSource for these specimens, and contact details of curatorial staff for specimens with restricted download (again, as shown on MorphoSource)

Dataset 6: Spreadsheet that contains age data, provenance data, and museum contacts for all fossil turtle specimens.

#####
TO PERFORM ANALYSES ON TURTLE LABYRINTH DATA

1. Store Datasets 2–5 & 7–8: these contain data to be read by the R-scripts
	- Dataset 2: specimen names, ID, measurements, ecological classifications, other info for turtle data
	- Dataset 3: Landmark data for 168 turtle specimens as csv files (3D coordinate data)
	- Dataset 4: information about which landmarks are used as sliding semilandmarks
	- Dataset 5: information about landmark colours used in one R script for deformation plots
	- Dataset 7: phylogenetic tree based on cal3-calibration, used for main analyses
	- Dataset 8: phylogenetic tree based on mil-calibration, used for sensitivity analyses
2. Open Datasets 9–13 & 16 in R to run different types of analyses. These scripts contain commands to read Datasets S1–6.
	- Dataset 9: performs GPA of membranous and endossoeus labyrinths, the does 2B-PLS to test landmarking scheme
	- Dataset 10: performs GPA & PCA, produces plots for Figure 1 of main text
	- Dataset 11: performs GPA & labyrinth shape regression analyses, produces Table 1 of main text
	- Dataset 12: performs GPA & size-corrected PCA, size-corrected shape regression analyses, produces Figure 2 of main text
	- Dataset 13: performs GPA & labyrinth size regression analyses, produces Table 2 of main text	
		- Dataset 14: output of 13, represents full version of table printed as Table 2 of main text
		- Dataset 15: as 14, but using different phylogenetic tree as a sensitivity analysis
	- Dataset 16: prints regression and residual plot, and performs ancestral state reconstructions; produces Figure 3 of main text
3. Note that Datasets 14–15 are an output list of the script Dataset 13


#####
TO PERFORM ANALYSES ON AMNIOTE LABYRINTH DATA

1. Store Datasets 17–19: these contain data to be read by the R-scripts
	- Dataset 17: specimen names, ID, measurements, ecological classifications, other info for amniote data
	- Dataset 18: Landmark data for 200 amniote specimens as csv files (3D coordinate data)
	- Dataset 19: information about which landmarks are used as sliding semilandmarks
2. Open Dataset 20 in R to run the amniote analysis. This script contains commands to read Datasets 17–19. Produces Figure 4 of main text.

#####
Dataset 21: This is the Source Data file for main text figures and supplementary figures that use numeric values plotted in graphs. The spreadsheet has several sheets, all for a figure or figure panel. These can be used to reproduce the figures directly (i.e., without having to run the R scripts, which produce the same figures).
