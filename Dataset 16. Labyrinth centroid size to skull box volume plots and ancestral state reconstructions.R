# Load packages

require( geomorph )
require( ape )
require( qpcR )
require( nlme )


	##Clear workspace
		rm( list = ls() )

	#Set working directory to where turtle landmark data (from supplemental folder: Dataset 3. Turtle labyrinth data) is stored
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

### The following reproduces the best model from the model comparison script

	GPA.data <- landmark.data.temp				
		tree.names <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "Tree_names" ] )
		fossil.temp <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "Specimen_type" ] )
			dimnames( GPA.data )[[ 3 ]][ !is.na( tree.names ) ] <- tree.names[ !is.na( tree.names ) ]
				duplicate.specimens <- which( is.na(tree.names) == TRUE )			
				fossil.specimens <- which( fossil.temp == "fossil")
					delete.these <- unique(c(duplicate.specimens, fossil.specimens))
						GPA.data <- GPA.data[,, - delete.these]
								
			GPA.labyrinth.extant <- gpagen( GPA.data , curves = sliders , ProcD = F )
				labyrinth.Csize <- GPA.labyrinth.extant$Csize
					labyrinth.Csize[ labyrinth.Csize > 5000 ] <- labyrinth.Csize[ labyrinth.Csize > 5000 ] / 1000
							#exclude because it has no complete cranial measurements
							labyrinth.Csize <- labyrinth.Csize[ - which(names(labyrinth.Csize) == "Trachemys_scripta") ]

		#Prepare tree that has same tips as the shape data blocks
			tree.temp <- drop.tip( tree , tree$tip.label[ ! tree$tip.label %in% names( labyrinth.Csize ) ] )

		#Prepare metadata that matches the labyrinth size data
			data.temp <- specimen.info[ specimen.info$Tree_names %in% tree.temp$tip.label , ]
				rownames( data.temp ) <-  data.temp$Tree_names
					data.temp <- data.temp[ tree.temp$tip.label , ]		

		#Extract vecological variable needed for best model. 
		
			notwebbed <- data.temp$Forelimb == "not webbed - 0"
				names( notwebbed ) <- rownames( data.temp )
			
		#Extract size variable variable needed for best model.		
				
			skull_box <- data.temp[ tree.temp$tip.label , "logV_mm3" ]
					names(skull_box) <- rownames( data.temp )

		#Make a data frame for analysss
		
			df <- data.frame( notwebbed, skull_box = skull_box, labyrinth_Csize = log10( labyrinth.Csize )[ tree.temp$tip.label ] , species = names(notwebbed) )

		#Reproduced best model 

			pgls.fit.best.model <- gls( labyrinth_Csize ~ skull_box + notwebbed , correlation= corPagel( 0.9, phy= tree.temp ), data = df )
			
				#warning message could be avoided by: pgls.fit.best.model <- gls( labyrinth_Csize ~ skull_box + notwebbed , correlation= corPagel( 0.9, phy= tree.temp, form = ~species ), data = df )

										
### The following uses the coefficients of the best model with the full dataset including fossils to get residual estimates		
### We also need the centroid sizes of all taxa to be included, thus we need to run another GPA with a larger dataset

	#Do GPA of labyrinth shape for all taxa available
			GPA.data <- landmark.data.temp
			
			#exclude all membranous specimens from main dataset
			
			membranous.specimens <- c(which(grepl("membranous", dimnames( GPA.data )[[ 3 ]] ) == TRUE))
				membranous.specimens.check <- dimnames( GPA.data )[[ 3 ]][membranous.specimens]
					GPA.data <- GPA.data[,, - membranous.specimens]	
				
				tree.names <- as.character( specimen.info[ dimnames( GPA.data )[[ 3 ]] , "Tree_names" ] )
					dimnames( GPA.data )[[ 3 ]][ !is.na( tree.names ) ] <- tree.names[ !is.na( tree.names ) ]			
			
			GPA.labyrinth.all <- gpagen( GPA.data , curves = sliders , ProcD = F )
			
				labyrinth.Csize.all <- GPA.labyrinth.all$Csize
					labyrinth.Csize.all[ labyrinth.Csize.all > 5000 ] <- labyrinth.Csize.all[ labyrinth.Csize.all > 5000 ] / 1000
							length(labyrinth.Csize.all)

		#DO PCA		
			PCA.labyrinth <- plotTangentSpace( GPA.labyrinth.all$coords , warpgrids = F )	
		
		
	#Prepare datasets in which tree, skull box volumes and labyrinth sizes pertain to the same set of taxa
			skull_box_full.temp <- specimen.info$logV_mm3
				names( skull_box_full.temp ) <- specimen.info[ , "Tree_names" ]
					skull_box_full.temp <- skull_box_full.temp[ !is.na( names( skull_box_full.temp ) ) ]
					skull_box_full.temp <- skull_box_full.temp[ !is.na( skull_box_full.temp ) ]

		#Prepares tree for ancestral state reconstructions
			tree.full.temp <- drop.tip( tree , tree$tip.label[ !tree$tip.label %in% names( labyrinth.Csize.all ) ] )
			tree.full.temp <- drop.tip( tree.full.temp , tree$tip.label[ !tree$tip.label %in% names( skull_box_full.temp ) ] )				
				tree.full.temp$edge.length[ which(tree.full.temp$edge.length == 0) ] <- 0.001

		#This is the labyrinth size and skull box data for all species included in the tree
			labyrinth_size_all <- labyrinth.Csize.all[ tree.full.temp$tip.label ]
			skull_box_full.temp <- skull_box_full.temp[ tree.full.temp$tip.label ]

	#Predict residuals for all species based on PGLS regression results
		predicted_labyrinth.CSize_full <- pgls.fit.best.model$coef[1] + pgls.fit.best.model$coef[2] * skull_box_full.temp
		residual_labyrinth.CSize_full <- log10(labyrinth_size_all) - predicted_labyrinth.CSize_full
		
		
	#Scale residuals from 0 to 1
		residuals.relative <- (as.numeric( residual_labyrinth.CSize_full ) - min( as.numeric( residual_labyrinth.CSize_full ))) / max(as.numeric( residual_labyrinth.CSize_full ) - min( as.numeric( residual_labyrinth.CSize_full )))
			names(residuals.relative) <- names(residual_labyrinth.CSize_full)

	
###		##Makes colour scale according to realtive labyrinth size
				colorScale <- colorRamp( c( "blue4", "blue" , "lightblue", "lightpink2", "red3", "deeppink4" ) )
					labyrinth_relative_size <- residuals.relative
						names( labyrinth_relative_size ) <- names( residuals.relative )
					labyrinth_relative_size <- c( labyrinth_relative_size , ace( labyrinth_relative_size , tree.full.temp )$ace )					
			BGs <- colorScale( labyrinth_relative_size )
				BG <- c()
					for( i in 1:nrow( BGs ) ) { BG[ i ] <- rgb( BGs[i,1] ,BGs[i,2] , BGs[i,3] , maxColorValue = 255 ) }
						names( BG ) <- names( labyrinth_relative_size )		


	# This makes a colour object in which specimens are coloured according to ecology
	
	plotting.habitat <- specimen.info[ c( match(tree.full.temp$tip.label, specimen.info$Species_ID) ) ,]$Plotting_habitat

	ecology.bg <- c( "dodgerblue3" , "lightblue3" , "tan" )
			names( ecology.bg ) <- c( "marine" , "freshwater" , "terrestrial" )
				ecology.bg <- ecology.bg[ as.character( plotting.habitat ) ]
					names(ecology.bg) <- tree.full.temp $tip.label
		
		##Makes tree plot
			
			dev.new(width=5, height=8)
			pdf("Labyrinth_Size_residuals_phylogeny.pdf", useDingbats=FALSE, width=5, height=10)
			plot.phylo( tree.full.temp , edge.col = BG[ tree.full.temp $edge[,2] ] , edge.width = 2.5 , cex = 0.4 , label.offset = 2 , no.margin = T, show.tip.label=F )
				#I'm putting the graph up a bit so the time scale fits underneath		
					par(xpd=T, mar= par()$mar + c(2,0,0,0))
						plot.phylo( tree.full.temp , edge.col = BG[ tree.full.temp$edge[,2] ] , edge.width = 2.5 , cex = 0.4 , label.offset = 2 , no.margin = F, show.tip.label=F )		
							
							#for numbers as tips: I am using the tip.label numbers of the TOTAL SUPERTREE WITH ALL TAXA, so that the taxon numbers
							#are consistent throughout different plots based on slightly differently pruned tree versions							
							numbers <- specimen.info[ match( tree.full.temp$tip.label ,  specimen.info$Tree_names ) , "Tree_node_number"  ]							
								tiplabels(numbers, frame="none", offset = 5.5, cex=0.4, font=2, bg = ecology.bg[ tree.full.temp $tip.label ], col = "black")						
							#for a time scale
							axisPhylo(side=1, cex.axis=0.5, pos=-1, lwd=2, mgp=c(3,0.3,0))
							#
							nodelabels( pch = 21 , bg = BG[ (length( tree.full.temp$tip.label ) + 1):(length( tree.full.temp$tip.label ) + tree.full.temp$Nnode ) ] , col = BG[ (length( tree.full.temp$tip.label ) + 1):(length( tree.full.temp$tip.label ) + tree.full.temp$Nnode ) ] )					
							tiplabels( pch = 21 , bg = BG[ tree.full.temp $tip.label ] , col = BG[ tree.full.temp $tip.label ] )
							#for including a clr legend
							require(plotrix)
							sortedBG <- colorScale(sort(labyrinth_relative_size[1:141]))
								BG2 <- c()
									for( i in 1:nrow( sortedBG ) ) { BG2[ i ] <- rgb( sortedBG[i,1] , sortedBG[i,2] , sortedBG[i,3] , maxColorValue = 255 ) }
										names( BG2 ) <- names( labyrinth_relative_size[1:141] )
							color.legend(90,18,150,21,legend = c("min", "max"), rect.col = BG2, cex=0.8, align="rb")
				dev.off()		
	
	#Csizes and box volumes for all specimens as objects that can be used for plotting	
		centroid.sizes.combined <- as.matrix(labyrinth_size_all)
		box.volumes.combined <- as.matrix(skull_box_full.temp)
	
	##This only plots the data used for tree plot:
	#	dev.new()
	#	plot( y = labyrinth_size_all , x = skull_box_full.temp , bty = "l" , log = "y" , xlab = "Skull box volume (mm^3)" , ylab = "Labyrinth centroid size (mm)" , pch = 21 , bg = BG , col = BG , cex = 2 )
	#	points(y = labyrinth_size_all , x = skull_box_full.temp , bty = "l" , log = "xy" , xlab = "Skull box volume (mm^3)" , ylab = "Labyrinth centroid size (mm)" , pch = 21 , bg = BG , col = BG , cex = 2)
	#	text( y = labyrinth_size_all , x = skull_box_full.temp , names( labyrinth_size_all ) , pos = 4 , cex = 0.6 )
	
		
	##The following scrip prepares the plot of all specimens for which data are available:
		
		#Extract Csizes and box volumes for additional specimens which are duplicates in terms of species representation			
		skull_box_extra <- specimen.info$logV_mm3
			names( skull_box_extra ) <- specimen.info[ , "Tree_names" ]
				skull_box_extra <- skull_box_extra[ !is.na( skull_box_extra ) ]	
				skull_box_extra <- skull_box_extra[ is.na( names( skull_box_extra ) ) ]
					extra_specimens <- c(which(is.na(specimen.info[ , "Tree_names" ])))
						drop <- which(is.na(specimen.info[extra_specimens, "logV_mm3"]))
						extra_specimens <- extra_specimens[-drop]
				names(skull_box_extra) <- specimen.info[ , "Specimen_name" ][extra_specimens]		
				skull_box_extra <- skull_box_extra[ !grepl("membranous", names(skull_box_extra))]
				skull_box_extra <- skull_box_extra[ !grepl("membraneous", names(skull_box_extra))]				
			
			#CSizes for additional specimens					
			labyrinth.CSize_extra <- labyrinth.Csize.all[names(skull_box_extra)]	
				
			## These could be added directly onto the plot of data used for PGLS:
				#points (skull_box_extra, labyrinth.CSize_extra, pch = 21, bg = "green", col="black", cex=2)		
					#text( y = labyrinth.CSize_extra , x = skull_box_extra , names(labyrinth.CSize_extra)  , pos = 4 , cex = 0.6 )	
								
		#Predict their residuals based on PGLS regression results
		predicted_labyrinth.CSize_extra <- pgls.fit.best.model$coef[1] + pgls.fit.best.model$coef[2] * skull_box_extra
		residual_labyrinth.CSize_extra <- log10(labyrinth.CSize_extra) - predicted_labyrinth.CSize_extra
			
		#Combine 'true' residuals and predicted residuals into a single object
		
		combined <- c(residual_labyrinth.CSize_full, residual_labyrinth.CSize_extra)
			
		#Scale combined residuals fro 0 to 1
		combined.residuals.relative <- (as.numeric( combined ) - min( as.numeric( combined ))) / max(as.numeric( combined ) - min( as.numeric( combined )))
				
			#Define color scheme
			BGs <- colorScale (combined.residuals.relative)
				BG <- c()
					for (i in 1:nrow( BGs )) {
						if (!is.na(BGs[i,1])) {
							BG[i] <- rgb( BGs[i,1], BGs[i,2], BGs[i,3], maxColorValue = 255 )
						} }
			
				names(BG) <- names(combined.residuals.relative)
				BG[is.na(BG)] <- "yellow" #colour is easy to catch in case an error happened
			
		#Combine Csizes and box volumes for all specimens into objects that can be used for plotting	
		centroid.sizes.combined2 <- rbind(as.matrix(centroid.sizes.combined), as.matrix(labyrinth.CSize_extra))	
		box.volumes.combined2 <- rbind(as.matrix(box.volumes.combined), as.matrix(skull_box_extra))
			
			#These could be projetced onto the existing plot, overrplotting the points with the correct color scheme
			#points (box.volumes.combined, centroid.sizes.combined, pch = 21, bg = BG, col="black", cex=2)
		
	# The following makes a colour element in which specimens are coloured according to ecology
		
	plotting.habitat.extra <- specimen.info[ c( match(names(labyrinth.CSize_extra), specimen.info$Specimen_name) ) ,]$Plotting_habitat

	ecology.bg.extra <- c( "dodgerblue3" , "lightblue3" , "tan" )
			names( ecology.bg.extra ) <- c( "marine" , "freshwater" , "terrestrial" )
				ecology.bg.extra <- ecology.bg.extra[ as.character( plotting.habitat.extra ) ]
					names(ecology.bg.extra) <- names(labyrinth.CSize_extra)
						ecology.bg <- c( ecology.bg , ecology.bg.extra )
			
		## Plot all together		
		dev.new()
		pdf("Centroid_size_box_volume_PGLS_ecology_labelled.pdf", useDingbats=FALSE)
		#plot( y = centroid.sizes.combined2 , x = box.volumes.combined2 , bty = "l" , log = "y" , xlab = "Skull box volume (mm^3)" , ylab = "Labyrinth centroid size (mm)" , pch = 21 , bg = BG , col = "black" , cex = 2 ) #alternative plot according to residual size
		plot( y = centroid.sizes.combined2 , x = box.volumes.combined2 , bty = "l" , log = "y" , xlab = "Skull box volume (mm^3)" , ylab = "Labyrinth centroid size (mm)" , pch = 21 , bg = ecology.bg , col = "black" , cex = 2 ) #alternative plot according to ecology
			abline( a = pgls.fit.best.model$coef[1], b = pgls.fit.best.model$coef[2] , lty = 2 , lwd = 1.5 , col = "grey" )	
			
				# Use full labels for points				
				#text( y = centroid.sizes.combined , x = box.volumes.combined , rownames(centroid.sizes.combined)  , pos = 4 , cex = 0.2 )
			
				# Text as numbers that equal number of tiplabel of full supertree
								
				text( y = labyrinth_size_all , x = skull_box_full.temp , numbers , cex = 0.4, font=2, offset=0.4 )
			
				#SWE prepare tiplabel objects for second batch of specimens. Numbers must match the numbers of the tree for the species for which each new datapoint represents a duplicate
					
					taxon.numbers <- specimen.info[ match( names(skull_box_extra) ,  specimen.info$Specimen_name ) , "Tree_node_number"  ]
								
				# Text as numbers for additional specimens; match numbers of corresponding speices from above 
				text( x = skull_box_extra , y = labyrinth.CSize_extra , taxon.numbers , cex = 0.4, font=2, offset=0.3 )
		
		dev.off()		
		
							
	###Residual plot
			
			groups <- as.character( specimen.info[ match( names( labyrinth_size_all ) , specimen.info$Tree_names ) , "Plotting_group" ] )
			
				#Combine some for easier plotting
					kinos <- which(groups == "Kinosternoidea")
						groups[ kinos ] <- "Chelydroidea"
					chelydrids <- c(which(groups == "Chelydridae"))
						groups[ chelydrids ] <- "Chelydroidea"
					xinj <- which(groups == "Xinjiangchelyidae")
						groups[ xinj ] <- "'Eucryptodires'"
					macros <- which(groups == "Macrobaenidae/Sinemydidae")
						groups[ macros ] <- "'Eucryptodires'"
					angolas <- which(groups == "Angolachelonia")
						groups[ angolas ] <- "'Eucryptodires'"
				
				#write data object for residual data containig family names
					data.part.one <- cbind(residual_labyrinth.CSize_full, groups)
					
				#define families for additional (non pgls) data 
					groups.part.two <- as.character( specimen.info[ match( names( residual_labyrinth.CSize_extra ) , specimen.info$Specimen_name ) , "Plotting_group" ] )
						kinos <- which(groups.part.two == "Kinosternoidea")
							groups.part.two[ kinos ] <- "Chelydroidea"
						chelydrids <- which(groups.part.two == "Kinosternoidea")
							groups.part.two[ chelydrids ] <- "Chelydroidea"
					
				#write data object for residual data containig family names
					data.parts.two <- cbind(residual_labyrinth.CSize_extra, groups.part.two)
					
				#combine datasets
					dataset.residual.plot <- rbind(data.part.one, data.parts.two)
							
						plot.names <- as.character(dataset.residual.plot[,2])
						res.values <- as.numeric(dataset.residual.plot[,1])
					
					#write is as dataframe		
					dataset.residual.plot <- data.frame(plot.names, res.values, rownames=rownames(dataset.residual.plot))

				#writes an object of numbers which correspond to taxon labels, matching tree tip labels
					number.labels <- c(numbers, taxon.numbers)	
					
			# Load package that allowed staggered axis labels with staxlab()
			require(plotrix)

		# Plot residual plot
			dev.new()
			pdf("pgls_CSize_residuals_numbered_ecology.pdf", useDingbats=FALSE)
				#Sorts clade categories in sequence that makes sense
				dataset.residual.plot$plot.names <- factor(dataset.residual.plot$plot.names, levels = c("early-diverging", "Perichelydia", "'Eucryptodires'", "Chelidae", "Pelomedusoides", "Trionychia", "Emysternia", "Geoemydidae", "Testudinidae",  "Chelydroidea", "Chelonioidea"))
					#Produces white boxplot as frame for later plotting arguments
					boxplot(dataset.residual.plot$res.values ~ dataset.residual.plot$plot.names, border="white", cex.axis=0.6, ylab="PGLS residuals", xlab="", names=F)
						#Plots ablines and points, using previous colour scale
						abline(h=0, col="grey")
						#points(dataset.residual.plot$res.values ~ dataset.residual.plot$plot.names, bg=BG, col="black", pch=21, cex=2) #plot according to residual size
						points(dataset.residual.plot$res.values ~ dataset.residual.plot$plot.names, bg=ecology.bg, col="black", pch=21, cex=2) #plot according to ecology
						#Includes taxon labels:
						#text( y = dataset.residual.plot$res.values , x = dataset.residual.plot$plot.names , dataset.residual.plot$rownames  , pos = 4 , cex = 0.2, font = 2 )
						#Includes taxon labels as numbers:
						text( y = dataset.residual.plot$res.values , x = dataset.residual.plot$plot.names , number.labels , cex = 0.4, font = 2 )
			staxlab(1,1:11,c("early stem turtles", "'Perichelydia'", "'Eucryptodires'", "Chelidae", "Pelomedusoides", "Trionychia", "Emysternia", "Geoemydidae", "Testudinidae",  "Chelydroidea", "Chelonioidea"), cex=0.7) 
			dev.off()					
					
### Make source file

	tree.full.temp	#is tree file used for tree
	tree.full.temp$tip.label # are tiplabels
	numbers	#are respective numbers used for each tip 
	residuals.relative #are relative residual sizes used for the ace function
	
	source.regression <- cbind(centroid.sizes.combined2, box.volumes.combined2)
		colnames(source.regression) <- c("log_labyrinth_centroid_sizes", "log_skull_box_volumes")
			residuals <- dataset.residual.plot$res.values
				source.file.3 <- cbind(source.regression, residuals)

	write.table(source.file.3, file = "Source_Data_file_3.txt")	
	