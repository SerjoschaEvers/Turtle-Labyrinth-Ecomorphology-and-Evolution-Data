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
			setwd( "INSERT DIRECTORY PATH" )
			sliders <- read.csv("Dataset 4. sliders.turtles.csv", row.names=1)
			colours <- as.character( read.csv("Dataset 5. landmark_colours.csv", row.names=1)[,1] )

		#Load specimen information
			setwd( "INSERT DIRECTORY PATH" )
			specimen.info <- read.csv( "Dataset 2. Specimen info.csv", header = TRUE )
				rownames( specimen.info ) <- specimen.info[ , "Specimen_name" ]

		#Load tree
			setwd( "INSERT DIRECTORY PATH" )		
			tree <- read.nexus( "Dataset 7. cal3tree.calibrated.txt" )
			alternative.tree <- read.nexus("Dataset 8. mbltree.calibrated.txt") 

		#Write GPA data objects for the membranous and endosseous data separately
			GPA.data <- landmark.data.temp
				membranous.specimens <- which(grepl("membranous", dimnames( GPA.data )[[ 3 ]]) == TRUE)
				endosseous.specimens <- which(grepl("endosseous", dimnames( GPA.data )[[ 3 ]]) == TRUE)
					GPA.data.membranous <- GPA.data[,, membranous.specimens ]
					GPA.data.endosseous <- GPA.data[,, endosseous.specimens ]
				
				#Delete inner loop landmarks
				ASC.loop.landmarks <- which(grepl("loop", dimnames( GPA.data.membranous )[[ 1 ]]) == TRUE)
				GPA.data.membranous <- GPA.data.membranous[-ASC.loop.landmarks,,]
				GPA.data.endosseous <- GPA.data.endosseous[-ASC.loop.landmarks,,]
					sliders <- sliders[- which( grepl("loop", rownames(sliders)) == TRUE) , ] 
						rows.to.modify <- c( which( grepl("LSC", rownames(sliders)) == TRUE) ,  which(grepl("PSC", rownames(sliders)) == TRUE ) )
						 sliders[rows.to.modify, ] <-sliders[rows.to.modify, ] -39		
			
			#Do GPA of both datasets			
			GPA.membranous <- gpagen( GPA.data.membranous , curves = sliders , ProcD = F )
			GPA.endosseous <- gpagen( GPA.data.endosseous , curves = sliders , ProcD = F )
			
			#DO PCA		###not necessary here, but in case I want to examine individual PCA plots
			PCA.membranous <- plotTangentSpace( GPA.membranous$coords , warpgrids = F )	
			PCA.endosseous <- plotTangentSpace( GPA.endosseous$coords , warpgrids = F )	
			
			# Change the names in the GPA object so they match tree names
			dimnames(GPA.membranous$coords)[[3]] <- as.character( specimen.info[ membranous.specimens , "Species_ID" ] )
			dimnames(GPA.endosseous$coords)[[3]] <- as.character( specimen.info[ endosseous.specimens , "Species_ID" ] )
				
			#Prepare tree that has same tips as the shape data blocks
			tree.temp <- drop.tip( tree , tree$tip.label[ ! tree$tip.label %in% dimnames( GPA.membranous$coords )[[ 3 ]] ] )
				plot(tree.temp)	
			
			##Run phylogenetic two-block partial least squares comparison of skull shape and labyrinth shape
				tbpls_fit <- phylo.integration( GPA.membranous$coords[ ,, tree.temp$tip.label ] , GPA.endosseous$coords[ ,, tree.temp$tip.label ] , phy = tree.temp )
			
			##You can look at the relative importance so the 2B-PLS axes here. You'll see that only the first few are important
					#These are proportions of 1.0
						round( (tbpls_fit$svd$d^2) / sum((tbpls_fit$svd$d^2)) , 3 )	#returns variances explained by all 2BPLS axes
						round( (tbpls_fit$svd$d^2) / sum((tbpls_fit$svd$d^2)) , 3 )[1:10]	#returns variances explained by first ten 2BPLS axes

			i = 1	#Use to select a 2BPLS axis for scatterplot and 3D shape deformations. "i=1" will select the first (most important) axis
				plot( tbpls_fit$XScores[,i] , tbpls_fit$YScores[,i] , xlim = range( tbpls_fit$XScores[,i] ) + c( 0 , 0.005 ), bty = "l" )
					text( tbpls_fit$XScores[,i] , tbpls_fit$YScores[,i] , labels = dimnames( tbpls_fit$XScores )[[ 1 ]] , cex = 0.6 , pos = 4 )
		