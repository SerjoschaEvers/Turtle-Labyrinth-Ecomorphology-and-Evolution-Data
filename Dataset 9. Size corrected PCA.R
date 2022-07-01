### Load packages
require( geomorph )
require( ape )

		#Clear workspace
		rm( list = ls() )

		#Set working directory to where turtle landmark data (from supplemental folder) is stored
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
			setwd( "~/Dropbox/Turtle DPhil Project MAC/Projects I - ongoing/Inner Ears Turtles/1 - Texts/5 - Nat Comms submission/1 - postreview/Github/v3.0" )
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

		#Do GPA of labyrinth shape for all taxa available
			GPA.data <- landmark.data.temp
				tree.names <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "Tree_names" ] )
				skull.box.temp <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "logV_mm3" ] )
					dimnames( GPA.data )[[ 3 ]][ !is.na( tree.names ) ] <- tree.names[ !is.na( tree.names ) ]	
						duplicate.specimens <- which( is.na(tree.names) == TRUE )
						no.skull.box <- which (is.na(skull.box.temp) == TRUE )
							delete.these <- unique(c(duplicate.specimens, no.skull.box))
								GPA.data <- GPA.data[,, - delete.these]
						
			GPA.labyrinth.all <- gpagen( GPA.data , curves = sliders , ProcD = F )
				labyrinth.Csize.all <- GPA.labyrinth.all$Csize
					labyrinth.Csize.all[ labyrinth.Csize.all > 5000 ] <- labyrinth.Csize.all[ labyrinth.Csize.all > 5000 ] / 1000

		#DO PCA		
			PCA.labyrinth <- plotTangentSpace( GPA.labyrinth.all$coords , warpgrids = F )	

		#Prepare tree that has same tips as the shape data blocks
			tree.temp <- drop.tip( tree , tree$tip.label[ ! tree$tip.label %in% names( GPA.labyrinth.all$Csize ) ] )
			
			#For tests with alternative calibration
			tree.temp.alternative <- drop.tip( alternative.tree , alternative.tree $tip.label[ ! alternative.tree $tip.label %in% names( GPA.labyrinth.all$Csize ) ] )

			#Examine trees
			plot(tree.temp, cex=0.4)
			plot(tree.temp.alternative, cex=0.4)

		##Make a version of the specimen data that matches the taxon sample
			data.temp <- specimen.info[ specimen.info$Tree_names %in% tree.temp$tip.label , ]
				rownames( data.temp ) <-  data.temp$Tree_names
					data.temp <- data.temp[ tree.temp$tip.label , ]
				
		##Size proxies
				
				skull_length.temp <- data.temp$Skull_length_mm
					names(skull_length.temp) <- rownames (data.temp) 	
						skull_length.temp <- skull_length.temp[ tree.temp$tip.label ]
					
									
				skull_width.temp <- data.temp$Skull_width_mm
					names(skull_width.temp) <- rownames (data.temp)
						skull_width.temp <- skull_width.temp[ tree.temp$tip.label ]
			
				
				skull_height.temp <- data.temp$Skull_height_mm
					names(skull_height.temp) <- rownames (data.temp)
						skull_height.temp <- skull_height.temp[  tree.temp$tip.label ]				
		
		##Skull geometry proxy
		
				skull_geometry.temp <- skull_height.temp / skull_width.temp
				
					#check frequency distribution
					hist(skull_geometry.temp)
				
					#check if these make sense
					which(skull_geometry.temp[] == max(skull_geometry.temp))
					which(skull_geometry.temp[] == min(skull_geometry.temp))
					
		
		#Make a data frame for analyses
				gdf <- geomorph.data.frame( shape = GPA.labyrinth.all$coords[ ,, tree.temp$tip.label ] ,  
						phy = tree.temp ,
						skull_box = data.temp[ tree.temp$tip.label , "logV_mm3" ] ,  						
						skull_geometry =  skull_geometry.temp)
												
		
		##In this script we're setting up all the combinations of explanatory variables for the right sizes of the models, 
		# and then running them all in a loop. This makes it easy to add regression models by extending the vector called "right.sides".
		# Here, we only use one model, which accounts for allometry and spatial constraints of braincase aspect ratio

		
			right.sides <- c(
			
			#following models test relations of size variables as correlates of shape, exploring allometric effects, plus the braincase aspect ratio
			"skull_box + skull_geometry" 
			
			)
			
				models <- paste( "shape ~" , right.sides )
					models <- lapply( models , as.formula )

		
		##Run Procrustes distance pGLS analyses (Adams 2014)
		
			procD.pgls.fits <- list()
				for( i in 1:length( models ) ) {
					procD.pgls.fits[[ i ]] <- procD.pgls( models[[ i ]] , phy = phy , SS.type = "II" , data = gdf )
				}
				
				##See summaries of procD.pgls results
					model.summaries <- lapply( procD.pgls.fits , summary )			
				
				#Print all coefficents to file
				capture.output(model.summaries, file = "Allometry_and_braincase_spatial_constraint_model.txt")


#After running the 'shape ~ skull size + braincase ratio' model:

	procD.model <- 	procD.pgls.fits[[ 1 ]]
	
	#Numbers for data point identification
		numbers <- specimen.info[ match( tree.temp $tip.label , specimen.info$Tree_names ) , "Tree_node_number"  ]	
	
	#Colour sheme for ecological categories
		ecology.bg <- c( "dodgerblue3" , "lightblue2" , "tan" )
			names( ecology.bg ) <- c( "marine" , "freshwater" , "terrestrial" )
				ecology.bg <- ecology.bg[ as.character( specimen.info[ , "Plotting_habitat" ][ match ( dimnames(procD.model$pgls.residuals)[[1]] , specimen.info$Species_ID ) ] ) ]

	#Panel A - Skull size scores
		dev.new(width = 8 , height = 4)
			close.screen( all.screens = T )
				split.screen( c( 1 , 2 ) )
		
		screen( 1 )	; par( mar = c( 4, 4, 1, 1 ) )
		size.scores <- plot ( procD.model, type="regression", reg.type="RegScore", predictor = procD.model$data$skull_box)$RegScore

			plot(size.scores ~ procD.model$data$skull_box, xlab="Skull box volume (mm^3)", ylab="Shape regression scores", pch=21, bg=ecology.bg, col = "darkgrey", lwd = 1.2, cex = 1.9)
				#text(size.scores ~ procD.model$data$skull_box, cex = 0.5 , labels = rownames(size.scores) )
				text(size.scores ~ procD.model$data$skull_box, cex = 0.3 , labels = numbers , col = "black" )
				legend("topleft", legend = c( "marine" , "freshwater" , "terrestrial" ) , pch = 21 , pt.bg = c( "dodgerblue3" , "lightblue2" , "tan" ) , cex = 1 , pt.cex = 1.5 , bty = "n" )

	#Panel B - Braincase ratio scores
		screen( 2 )	; par( mar = c( 4, 4, 1, 1 ) )
		bc.scores <- plot ( procD.model, type="regression", reg.type="RegScore", predictor=procD.model$data$skull_geometry)$RegScore

			plot(bc.scores ~ procD.model$data$skull_geometry, xlab="Braincase aspect ratio", ylab="Shape regression scores", pch=21, bg= ecology.bg, col = "darkgrey", lwd = 1.2, cex = 1.9)
			#text(bc.scores ~ procD.model$data$skull_geometry, cex = 0.5 , labels = rownames(size.scores) )
			text(bc.scores ~ procD.model$data$skull_geometry, cex = 0.3 , labels = numbers , col = "black" )

	#Panel C - PCA with procD residuals
		dev.new( width = 8 , height = 4 )
			close.screen( all.screens = T )
				split.screen( c( 1 , 2 ) )
		
		PCA.resid <- prcomp(procD.model$pgls.residuals)

		X=1
		Y=2

		screen( 1 )	; par( mar = c( 4, 4, 1, 1 ) )
			plot ( PCA.resid$x [,c(X,Y)] , bty="l", pch=21, bg = ecology.bg , col="darkgrey" , cex = 2 , 
				xlab = paste0("PC",X," (",round((PCA.resid$sdev^2 / sum(PCA.resid$sdev^2))*100,2)[X], "%",")") ,
				ylab = paste0("PC",Y," (",round((PCA.resid$sdev^2 / sum(PCA.resid$sdev^2))*100,2)[Y], "%",")")		)
				#text( PCA.resid$x [,c(X,Y)] , cex = 0.5 , labels = rownames(PCA.resid$x) )
				text(PCA.resid$x [,c(X,Y)] , cex = 0.3 , labels = numbers , col = "black" )
		
		Y=3
		
		screen( 2 )	; par( mar = c( 4, 4, 1, 1 ) )
			plot ( PCA.resid$x [,c(X,Y)] , bty="l", pch=21, bg = ecology.bg , col="darkgrey" , cex = 2 , 
				xlab = paste0("PC",X," (",round((PCA.resid$sdev^2 / sum(PCA.resid$sdev^2))*100,2)[X], "%",")") ,
				ylab = paste0("PC",Y," (",round((PCA.resid$sdev^2 / sum(PCA.resid$sdev^2))*100,2)[Y], "%",")")		)
				#text( PCA.resid$x [,c(X,Y)] , cex = 0.5 , labels = rownames(PCA.resid$x) )
				text(PCA.resid$x [,c(X,Y)] , cex = 0.3 , labels = numbers , col = "black" )

##Deformation plots
	
	#Plot shape deformation for "skull size" and "braincase ratio" variables
		#Get coefficients from 'procD.pgls' object
		coefficients.temp <- procD.model$pgls.coefficients 

	#Optional: Get 1st and 3rd quartiles of continuous variables to avoid potential outliers (minimum and maximum values) -> exchange c(0,1) with c(0.25,0.75)
	#Currently it takes the minimum and maximum values without averaging over a quartile

		skull.size.quart <- quantile(procD.model$data$skull_box,prob=c(0,1))
		bc.ratio.quart <- quantile(procD.model$data$skull_geometry,prob=c(0,1))
	
	##Run these next two blocks separately for deformation plots for skull size and braincase aspect ratio
		#Allometry variation (skull size)
			variables1 <- list( Intercept = 1  , skull_box = min(skull.size.quart) , skull_geometry = mean(procD.model$data$skull_geometry) )
			variables2 <- list( Intercept = 1  , skull_box = max(skull.size.quart) , skull_geometry = mean(procD.model$data$skull_geometry) )

		#Braincase ratio variation
			variables1 <- list( Intercept = 1  , skull_box = mean(procD.model$data$skull_box) , skull_geometry = min(bc.ratio.quart) )
			variables2 <- list( Intercept = 1  , skull_box = mean(procD.model$data$skull_box) , skull_geometry = max(bc.ratio.quart) )
	
	##All following lines must be run for each block separatly: 
	
		# Change the names of these lists to match the names in your 'coefficients.temp' object

			names( variables1 ) <- rownames ( coefficients.temp )
			names( variables2 ) <- rownames ( coefficients.temp )

		# Create lists to which you will save the values calculated for each landmark in each shape deformation	
			coefficient.list.temp1 <- list()
			coefficient.list.temp2 <- list()

		# Calculate predicted shapes (min and max)

			# (step 1)
				for( row.temp in 1:nrow( coefficients.temp ) ) {
 					 coefficient.list.temp1[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables1[[ rownames( coefficients.temp )[ row.temp ] ]]
 					 coefficient.list.temp2[[ row.temp ]] <- coefficients.temp[ row.temp , ] * variables2[[ rownames( coefficients.temp )[ row.temp ] ]]
					}

			# (step 2) 
				shape1 <- matrix( apply( do.call( rbind , coefficient.list.temp1 ) , 2 , sum ) , ncol = 3 , byrow = T ) 
				shape2 <- matrix( apply( do.call( rbind , coefficient.list.temp2 ) , 2 , sum ) , ncol = 3 , byrow = T )

			# define the colours you want
				point.distance.scale <- colorRamp(c("lightgrey" ,'red')) 

			point.distances <- c()
				Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 } #Euclidean distance function

			for ( i in 1:nrow(shape2) ){
    
   				 # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
   				 point.distances[i] <- Edist(shape1[i,] , shape2[i,]) 
				}	

			#point.distances

			# normalise the distances so they range from 0 to 1
				point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

			# and then you're able to create a colorRamp that goes from grey (0) to red (1)
				point.colours <- point.distance.scale(point.distances.norm)

	# 3D plot of shape deformation

		sliders.temp <- as.matrix(sliders)

		open3d()
			plot3d(shape1, size = 10, col='lightgrey' , box = "n" , aspect="iso" )
				lines.plot.temp <- rbind( sliders.temp[ , 1:2 ] , sliders.temp[ , 2:3 ] )
					for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape1[ lines.plot.temp[j,] , ] , lwd = 2 , col = 'lightgrey' ) }

			points3d(shape2, size = 15, col=rgb(point.colours,maxColorValue = 255) )
				lines.plot.temp <- rbind( sliders.temp[ , 1:2 ] , sliders.temp[ , 2:3 ] )
					for( j in 1:nrow( lines.plot.temp ) ) { lines3d( shape2[ lines.plot.temp[j,] , ] , lwd = 3 , col = rgb(point.colours,maxColorValue = 255) ) }

	rgl.snapshot( "Skull size deformation lateral.png" )
	rgl.snapshot( "Skull size deformation dorsal.png" )	
	rgl.snapshot( "Braincase aspect ratio deformation lateral.png" )	
	rgl.snapshot( "Braincase aspect ratio deformation dorsal.png" )												 
		