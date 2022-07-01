### Load packages

	require( geomorph )
	require( stringdist )
	require( ape )
	require( nlme )

###Clear workspace

	rm( list = ls() )

###Set working directory where amniote data (Dataset 13. Amniote labyrinth data) is stored

	setwd( "INSERT DIRECTORY PATH" )
	
	## Load landmark data
		
			temp.file <- list.files(pattern = ".csv") 			landmark.data <- lapply (temp.file, read.csv, row.names=1)
				names(landmark.data) <- gsub(".csv","",temp.file)
		
			landmark.data.temp <- array(as.numeric(unlist(landmark.data)), dim = c(236, 3, 200)) #number of landmarks, number of dimensions (3D coordinates), number of specimens
				dimnames(landmark.data.temp)[[3]] <- gsub(".csv","",temp.file)
				dimnames(landmark.data.temp)[[1]] <- rownames(landmark.data[[1]])
				dimnames(landmark.data.temp)[[2]] <- c("x","y","z")
	
	## Load slider info
		
		setwd( "INSERT DIRECTORY PATH" )
		sliders <- read.csv("Dataset 15. sliders.amniotes.csv", row.names=1)

	#Load specimen information
		setwd( "INSERT DIRECTORY PATH" )
		specimen.info <- read.csv( "Dataset 13. specimen.info.amniotes.csv", header = TRUE )
			rownames( specimen.info ) <- specimen.info[ , "Project_name" ]


	
### Define skull lengths (specimen.info has a skull length column)	
	
	postrostral_length <- specimen.info$Postrostral_skull_length_mm
		names(postrostral_length) <- specimen.info$Species
			postrostral_length[ postrostral_length > 2000 ] <- postrostral_length[ postrostral_length > 2000 ] /1000				
		
### Prepare labyrinth CSize variable	

		#Do GPA of labyrinth shape for all taxa available
			GPA.data <- landmark.data.temp					
			
			GPA.labyrinth.all <- gpagen( GPA.data , curves = sliders , ProcD = F )
				labyrinth.Csize.all <- GPA.labyrinth.all$Csize
					labyrinth.Csize.all[ labyrinth.Csize.all > 3000 ] <- labyrinth.Csize.all[ labyrinth.Csize.all > 3000 ] / 1000
						names(labyrinth.Csize.all) <- specimen.info[ match(names(labyrinth.Csize.all), specimen.info$Project_name), "Species" ] 
							which(labyrinth.Csize.all == max(labyrinth.Csize.all))

###Taxonomic categories for plotting	

	plot.group <- as.character( specimen.info[, "Clade" ] )
		names(plot.group) <- specimen.info$Species
			unique(plot.group)
						
			plot.group[which(names(plot.group) == "Aldabrachelys_gigantea")] <- "Testudinidae"
			plot.group[which(names(plot.group) == "Gopherus_agassizii")] <- "Testudinidae"
			plot.group[which(names(plot.group) == "Gopherus_flavomarginatus")] <- "Testudinidae"
			plot.group[which(names(plot.group) == "Manouria_impressa")] <- "Testudinidae"
			plot.group[which(names(plot.group) == "Testudo_marginata")] <- "Testudinidae"
			plot.group[which(names(plot.group) == "Homopus_areolatus")] <- "Testudinidae"
			
	good.colours <- function( n ) {
		c(   "aquamarine2" ,  "darkorchid2" ,  "dodgerblue2", "darkred" , "darkolivegreen4", "orange2", "deeppink2"   )[ 1:n ]
		}
	
		#For colour assignment
			specimen.order <- c( "Mammalia" , "Aves" , "Crocodylia" , "Lepidosauria" , "Archosauriformes", "Testudinata" , "Testudinidae" )
	
			#Assign colours according to clades
				group.bg <- good.colours( length( unique( plot.group ) ) )[ match( plot.group , specimen.order ) ]
					names( group.bg ) <- plot.group

###Makes numbers for indentification of taxa in plot							

	taxon.numbers <- seq( 1:length(rownames(specimen.info)) )

###Writes dataframe for plotting

	df <- data.frame( labyrinth_Csize = log10( labyrinth.Csize.all ) , postrostral_length = log10( postrostral_length ) )
	
	#export raw data to submit as supplement
	df.export <- data.frame( labyrinth_Centroid_size =  labyrinth.Csize.all  , postrostral_skull_length =  postrostral_length  )
	#Print table to file	
		write.table(df.export, file = "Tetrapod_relative_eas_size_data.csv", sep=" ")
				
	dev.new()
		plot( x = postrostral_length , y = labyrinth.Csize.all , bty = "l" , log = "xy" , xlab = "Postrostral skull length (mm)" , ylab = "Labyrinth centroid size (mm)" , pch = 21 , bg = group.bg , col = "black" , cex = 2 )
		# Text as numbers, which match sequence of specimens in specimen.info 
				text( x = postrostral_length , y = labyrinth.Csize.all , taxon.numbers , cex = 0.4, font=2, offset=0.3 )
		legend( "topleft" , legend = specimen.order , pch = 21 , pt.bg = good.colours( length( unique( plot.group ) ) ) , cex = 1 , pt.cex = 2.5 , bty = "n" )

