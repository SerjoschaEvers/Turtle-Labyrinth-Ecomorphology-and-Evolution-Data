### Load packages

require( geomorph )
require( ape )

	#Clear workspace
		rm( list = ls() )

		#Set working directory to where turtle landmark data (from supplemental folder: Dataset 2. Turtle labyrinth data) is stored
			setwd( "INSERT DIRECTORY PATH" )

		#Load landmark data and transform to correct format for GPA commands
		
			temp.file <- list.files(pattern = ".csv")
		
			landmark.data <- lapply (temp.file, read.csv, row.names=1)
				names(landmark.data) <- gsub(".csv","",temp.file)
		
			landmark.data.temp <- array(as.numeric(unlist(landmark.data)), dim = c(123, 3, 184)) #number of landmarks, number of dimensions (3D coordinates), number of specimens
				dimnames(landmark.data.temp)[[3]] <- gsub(".csv","",temp.file)
				dimnames(landmark.data.temp)[[1]] <- rownames(landmark.data[[1]])
				dimnames(landmark.data.temp)[[2]] <- c("x","y","z")

		#Load slider information and colour information for deformation plots
			setwd( "INSERT DIRECTORY PATH" )
			sliders <- read.csv("Dataset 3. sliders.turtles.csv", row.names=1)
			colours <- as.character( read.csv("Dataset 4. landmark_colours.csv", row.names=1)[,1] )

		#Load specimen information
			setwd( "INSERT DIRECTORY PATH" )
			specimen.info <- read.csv( "Dataset 1. Specimen info.csv", header = TRUE )
				rownames( specimen.info ) <- specimen.info[ , "Specimen_name" ]

		#Load tree
			setwd( "INSERT DIRECTORY PATH" )		
			tree <- read.nexus( "Dataset 5. cal3tree.calibrated.txt" )
			alternative.tree <- read.nexus("Dataset 6. mbltree.calibrated.txt")

	##GPA analysis and PCA

		#Do GPA of labyrinth shape for all taxa available
			GPA.data <- landmark.data.temp
			
			#exclude all membranous specimens from main dataset
			
			membranous.specimens <- c(which(grepl("membranous", dimnames( GPA.data )[[ 3 ]] ) == TRUE))
				membranous.specimens.check <- dimnames( GPA.data )[[ 3 ]][membranous.specimens]
					GPA.data <- GPA.data[,, - membranous.specimens]	
				
				tree.names <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "Tree_names" ] )
					dimnames( GPA.data )[[ 3 ]][ !is.na( tree.names ) ] <- tree.names[ !is.na( tree.names ) ]			
			GPA.labyrinth.all <- gpagen( GPA.data , curves = sliders , ProcD = F )

		#DO PCA		
			PCA.labyrinth <- plotTangentSpace( GPA.labyrinth.all$coords , warpgrids = F )
				superposed.lms <- GPA.labyrinth.all		


#Extract order of species in PCA data
datarownames <- rownames(PCA.labyrinth$pc.scores)
	get.two <- function( X ) { X[1:2] }
		speciesPCAorder <- unlist( lapply( lapply( strsplit( datarownames , "_" ) , get.two ) , paste , collapse = "_" ) )

			
		#Make reduced info file with columns relevant for plotting	
		reduced.info <- specimen.info[, c("Species_ID", "Plotting_group", "Plotting_habitat", "Replicates_endosseous", "Tree_node_number")]			
			
			#Delete all entries of specimens based on membranous labyrinths
			reduced.info <- reduced.info[ - c(which(grepl("membranous", rownames( reduced.info )) == TRUE)) , ]
			
				#Check if these are in the same order/sequence
				speciesPCAorder == reduced.info$Species_ID
				#These should now contain the same number of species
				length(reduced.info[, "Species_ID"]) == length(speciesPCAorder)

 	
	#Extract families from reduced info sheet, which are now ordered according the PCA data 	
	plot.group <- as.character( reduced.info[, "Plotting_group" ] )
		names(plot.group) <- rownames(reduced.info)
			unique(plot.group)

	#Colour sheme from Roger
		good.colours <- function( n ) {
		c( "firebrick2" ,  "seagreen3"  , "aquamarine2" , "slategray3" , "darkorchid2" , "cadetblue3" , "darkred" , "coral2" , "chartreuse3" , "grey48" , "orange2" , "springgreen2" ,  "salmon2" , "dodgerblue2" )[ 1:n ]
		}
	
		#For colour assignment
			specimen.order <- c( "early-diverging" , "Perichelydia" , "Xinjiangchelyidae" , "Macrobaenidae/Sinemydidae" , "Angolachelonia" , "Chelidae" , "Pelomedusoides" , "Trionychia" , "Testudinidae" , "Geoemydidae", "Emysternia", "Chelydridae" , "Kinosternoidea" , "Chelonioidea" )
	
			#Assign colours according to families
				group.bg <- good.colours( length( unique( plot.group ) ) )[ match( plot.group , specimen.order ) ]
					names( group.bg ) <- plot.group
							
					#Tests if colours indeed match families (works)
					#reduced.info[,"Colours"]<-group.bg
					#reduced.info[,-c(3:4)]

	#Colour sheme for ecological categories
		ecology.bg <- c( "dodgerblue3" , "lightblue2" , "tan" )
			names( ecology.bg ) <- c( "marine" , "freshwater" , "terrestrial" )
				ecology.bg <- ecology.bg[ as.character( reduced.info[ , "Plotting_habitat" ] ) ]


	#extract numbers for point labels on plot
		numbers.for.plot <- reduced.info[, "Tree_node_number"]
					

##Plot PC space
	PCA.results <- PCA.labyrinth
		labelled = T	#Select if you want taxon names to appear on plot
		labelled = F	#Select if you do not want taxon names to appear on plot

	#Choose which PC axes to show and create axis labels that state the proportion of variance explained
		X <- "PC1"
			XLAB <- paste( X , " (" , round( PCA.results$pc.summary$importance[ 2 , X ] , 3 ) * 100 , "%" , ")" , sep = "")
		Y <- "PC2"
			YLAB <- paste( Y , " (" , round( PCA.results$pc.summary$importance[ 2 , Y ] , 3 ) * 100 , "%" , ")" , sep = "")
		Z <- "PC3"
			ZLAB <- paste( Z , " (" , round( PCA.results$pc.summary$importance[ 2 , Z ] , 3 ) * 100 , "%" , ")" , sep = "")			

	dev.new( width = 8 , height = 8 )
			close.screen( all.screens = T )
		split.screen( c( 2 , 2 ) )
			screen( 1 )	; par( mar = c( 4, 4, 1, 1 ) )
				#I could use plot( PCA.results$pc.scores[ , c(X,Y) ] , parameters), but wanted to flip the PC1 axis so that my previous descriptive text is still true				
				plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , bty = "n" , pch = 21 , bg = group.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = YLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Y] ) + c( 0 , 0.07 ) , col = "darkgrey" )				
					if(labelled) {text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , rownames( PCA.results$pc.scores ) , cex = 0.2 )}
					text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , labels = numbers.for.plot , cex = 0.3 )
						legend( "topleft" , legend = specimen.order , pch = 21 , pt.bg = good.colours( length( unique( plot.group ) ) ) , cex = 0.5 , pt.cex = 1.5 , bty = "n" )
			screen( 2 )	; par( mar = c( 4, 4, 1, 1 ) )
				plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , bty = "n" , pch = 21 , bg = group.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = ZLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Z] ) + c( 0 , 0.07 )  , col = "darkgrey" )
					if(labelled) {text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , rownames( PCA.results$pc.scores ) , cex = 0.2 )}
					text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , labels = numbers.for.plot , cex = 0.3 )
			screen( 3 )	; par( mar = c( 4, 4, 1, 1 ) )
			plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , bty = "n" , pch = 21 , bg = ecology.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = YLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Y] ) + c( 0 , 0.07 ) , col = "darkgrey" )
					if(labelled) {text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , rownames( PCA.results$pc.scores ) , cex = 0.2 )}
					text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , labels = numbers.for.plot , cex = 0.3 )
						legend( "topleft" , legend = c( "marine" , "freshwater" , "terrestrial" ) , pch = 21 , pt.bg = c( "dodgerblue3" , "lightblue2" , "tan" ) , cex = 0.5 , pt.cex = 1.5 , bty = "n" )
			screen( 4 )	; par( mar = c( 4, 4, 1, 1 ) )
				plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , bty = "n" , pch = 21 , bg = ecology.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = ZLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Z] ) + c( 0 , 0.07 )  , col = "darkgrey" )
					if(labelled) {text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , rownames( PCA.results$pc.scores ) , cex = 0.2 )}
					text( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , labels = numbers.for.plot , cex = 0.3 )


# Plot individual families against all other points

	#Define colours for single families and the legend
	i <- specimen.order[1]	

		single.group.bg <- group.bg
			single.group.bg[ which(names(single.group.bg[]) != i ) ] <- "grey"
	
		for.labels <- c(which(single.group.bg[] != "grey"))
	
		legend.single.group.text <- c(i, "other")
		legend.single.group.bg <- c(single.group.bg[[i]], "grey")
				
		for.points.x <-  -PCA.results$pc.scores[ , X ][which(plot.group == i) ]
		for.points.y <- PCA.results$pc.scores[ , Y ][which(plot.group == i) ]
		for.points.z <- PCA.results$pc.scores[ , Z ][which(plot.group == i) ]		

	#Plot
	dev.new( width = 8 , height = 4 )
				close.screen( all.screens = T )
			split.screen( c( 1 , 2 ) )
				screen( 1 )	; par( mar = c( 4, 4, 1, 1 ) )
					plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Y ] , bty = "n" , pch = 21 , bg = single.group.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = YLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Y] ) + c( 0 , 0.07 ) , col = "darkgrey" )						
					points( for.points.x , for.points.y , bty = "n" , pch = 21 , bg = group.bg[which(names(group.bg[]) == i )], lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = YLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Y] ) + c( 0 , 0.07 ) , col = "darkgrey" )							
						if(labelled) {text( for.points.x , for.points.y , rownames( PCA.results$pc.scores )[for.labels] , cex = 0.2 )}
						text( for.points.x , for.points.y , labels = numbers.for.plot[for.labels] , cex = 0.3 )
							legend( "topleft" , legend = legend.single.group.text , pch = 21 , pt.bg = legend.single.group.bg , cex = 0.5 , pt.cex = 1.5 , bty = "n" )
				screen( 2 )	; par( mar = c( 4, 4, 1, 1 ) )
					plot( -PCA.results$pc.scores[ , X ], PCA.results$pc.scores[ , Z ] , bty = "n" , pch = 21 , bg = single.group.bg, lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = ZLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Z] ) + c( 0 , 0.07 )  , col = "darkgrey" )					
					points( for.points.x , for.points.z  , bty = "n" , pch = 21 , bg = group.bg[which(names(group.bg[]) == i )], lwd = 1.2 , cex = 1.9 , xlab = XLAB , ylab = ZLAB , cex.lab = 0.9 , cex.axis = 0.8 , ylim = range( PCA.results$pc.scores[,Z] ) + c( 0 , 0.07 )  , col = "darkgrey" )					
						if(labelled) {text( for.points.x , for.points.z , rownames( PCA.results$pc.scores )[for.labels] , cex = 0.2 )}
						text( for.points.x , for.points.z , labels = numbers.for.plot[for.labels] , cex = 0.3 )




					
####Plot PC shapes
	##It all runs much faster if you omit the lines that plot lines connecting landmarks together

labelled = FALSE
		lines.plot <- rbind( as.matrix(sliders[,c(1,2)]) , as.matrix(sliders[,c(2,3)]) )

	open3d()
		par3d(windowRect = c(0,0,1400,900))
			Sys.sleep(1)
				mfrow3d(nr = 3, nc = 4, byrow = TRUE, sharedMouse = TRUE)
	#to match reversed PC1 axis I plit the max shape first and label it as min			
	choose <- "PC1max"	#choose <- "PC1max"; reversed
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main="PC1min")} )
				aspect3d( "iso" )							
			next3d()
	#to match reversed PC1 axis I plit the min shape second and label it as max	
	choose <- "PC1min"	#choose <- "PC1min"; reversed so this fits with my reversed axes in the plot
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main="PC1max")} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC4min"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC4max"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC2min"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC2max"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
	choose <- "PC5min"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC5max"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )				
			next3d()
	choose <- "PC3min"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC3max"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC6min"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )
			next3d()
	choose <- "PC6max"
		plot3d( PCA.labyrinth$pc.shapes[[choose]] , col = colours , size = 10 , box = "n")
			for( i in 1:nrow( lines.plot ) ) { lines3d( PCA.labyrinth$pc.shapes[[choose]][ lines.plot[i,] , ] , lwd = 2 , col = "darkgrey" ) }
				points3d( superposed.lms$consensus , col = "grey" , size = 5 , box = "n")
			if( labelled) { title3d( main = choose, line=1 , cex = 1) }
								bgplot3d( { plot.new() ;title(main=choose)} )
				aspect3d( "iso" )					
				
				
##Screenshots for PC shapes. 
				
	rgl.snapshot( "PC shapes 1-6 lateral view turtles.png" )
	rgl.snapshot( "PC shapes 1-6 dorsal view turtles.png" )
	rgl.snapshot( "PC shapes 1-6 anterior canal view turtles.png" )
	rgl.snapshot( "PC shapes 1-6 posterior canal view turtles.png" )			
	
#			