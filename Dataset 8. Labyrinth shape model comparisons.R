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


		##Make explanatory variables
		
			##Habitat ecology	

			marine.all <-  data.temp$Habitat_general == "marine" | data.temp$Fossil_marine == "yes"
				names( marine.all ) <- rownames( data.temp )	
					marine.all <- marine.all[ tree.temp$tip.label ]						
						marine.all[which(is.na(marine.all))] <- "FALSE"
							which(marine.all == TRUE)
			
			marine <-  data.temp$Habitat_general == "marine" 
				names( marine ) <- rownames( data.temp )	
					marine <- marine[ tree.temp$tip.label ]					
						marine[which(is.na(marine))] <- "FALSE"
							which(marine == TRUE)
				
			freshwater <-  data.temp$Plotting_habitat == "freshwater"  
				names( freshwater ) <- rownames( data.temp )	
					freshwater <- freshwater[ tree.temp$tip.label ]
						freshwater[which(is.na(freshwater))] <- "FALSE"
							which(freshwater == TRUE)
				
			terrestrial <-  data.temp$Plotting_habitat == "terrestrial" 
				names( terrestrial ) <- rownames( data.temp )
					terrestrial <- terrestrial[ tree.temp$tip.label ]
						terrestrial[which(is.na(terrestrial))] <- "FALSE"
							which(terrestrial == TRUE)
			
			##Neck categories
		
			no_plane <-  data.temp$Retraction_type == "none" 
				names( no_plane ) <- rownames( data.temp )
					no_plane <- no_plane[ tree.temp$tip.label ]
						no_plane[which(is.na(no_plane))] <- "FALSE"
							which(no_plane == TRUE)
			
			vertical <-  data.temp$Retraction_type == "vertical" 
				names( vertical ) <- rownames( data.temp )
					vertical <- vertical[ tree.temp$tip.label ]
						vertical[which(is.na(vertical))] <- "FALSE"
							which(vertical == TRUE)
							
			horizontal <-  data.temp$Retraction_type == "sideways" 
				names( horizontal ) <- rownames( data.temp )
					horizontal <- horizontal[ tree.temp$tip.label ]
						horizontal[which(is.na(horizontal))] <- "FALSE"
							which(horizontal == TRUE)
							
			incomplete_retr <-  data.temp$Retratction_ability == "incomplete" 
				names( incomplete_retr ) <- rownames( data.temp )
					incomplete_retr <- incomplete_retr[ tree.temp$tip.label ]
						incomplete_retr[which(is.na(incomplete_retr))] <- "FALSE"
							which(incomplete_retr == TRUE)
			
			full_retr <-  data.temp$Retratction_ability == "full" 
				names( full_retr ) <- rownames( data.temp )
					full_retr <- full_retr[ tree.temp$tip.label ]
						full_retr[which(is.na(full_retr))] <- "FALSE"
							which(full_retr == TRUE)
				
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
					
		
		#Make a big data frame for analyses
				gdf <- geomorph.data.frame( shape = GPA.labyrinth.all$coords[ ,, tree.temp$tip.label ] ,  
						phy = tree.temp ,
						marine = marine , freshwater = freshwater , terrestrial = terrestrial , marine.all = marine.all, 
						no_plane = no_plane, vertical = vertical, horizontal = horizontal,
						incomplete_retr = incomplete_retr, full_retr = full_retr,					 
						skull_length = log10( skull_length.temp ) , skull_width = log10(skull_width.temp) , skull_height = log10(skull_height.temp) , 
						skull_box = data.temp[ tree.temp$tip.label , "logV_mm3" ] ,  
						labyrinth_Csize = log10( labyrinth.Csize.all )[ tree.temp$tip.label ] ,
						skull_geometry =  skull_geometry.temp)
												
		
		##In this script we're setting up all the combinations of explanatory variables for the right sizes of the models, 
		# and then running them all in a loop. This makes it easy to add regression models by extending the vector called "right.sides".
		
		#as the model building process was iterative, several models that were initially explored are muted below. 
		#Models that are active are those reported in the table <shape_models_incl_fossils>
		
			right.sides <- c(
			
			#following models test relations of size variables as correlates of shape, exploring allometric effects
			"skull_length" , "skull_width",  "skull_height", "skull_box" , "labyrinth_Csize" , 
				# -> skull height performs best (R2), followed by skull box.
			
			#following model tests relations of skull geometry as correlates of shape 
			"skull_geometry" ,  
				# -> significant and explains much of the variance
			
			#following models test independent effect of skull size and labyrinth size
			"skull_length + labyrinth_Csize" , "skull_box + labyrinth_Csize" ,   
				#-> both models show significant indipendent effect, skull box performs better
			
			#following models test independent effect of skull size and skull geometry
			"skull_length + skull_geometry" , "skull_box + skull_geometry" , "skull_height + skull_geometry",	
				#-> all models significant. in the height+geometry model, proportion of variance explained is near equal between both variables
				#-> in other models, more variance is explanation is attributed to geometry than skull size
			
			#following models test independent effect of skull geometry and labyrinth size
			"skull_geometry + labyrinth_Csize" , 
				# -> also significant
			
			#following models test independent effects of skull size and skull geometry and labyrinth size
			"skull_length + skull_geometry + labyrinth_Csize" , #slightly worse in R2 than below model
			"skull_box + skull_geometry + labyrinth_Csize" , #slightly better in R2 than above model
				#-> all independent effects are important
						
			#following models test non-independent effects, i.e. hypothesis that taxa with prop. larger labyrinths in relation to skull size have different laby shapes
			"skull_length * labyrinth_Csize" , #
			"skull_box * labyrinth_Csize" ,  #interaction term significant
			
			#
			"skull_box * labyrinth_Csize + skull_geometry" ,
				#->interaction term remains significant 
			
			#following models test non-independent effects, i.e. hypothesis that taxa with higher/broader skulls in relation to skull size have different laby shapes
			"skull_length * skull_geometry" , #slightly worse in R2 than below model
			"skull_box * skull_geometry" , #slightly better in R2 than above model
				#-> interaction term is significant
								   
			#the following models ask: do terrestrial turtles have a different mean labyrinth shape than non-terrestrial trutles?
			# -> initial analyses show:  only freshwater and terretrial are relevant
			"terrestrial" , "freshwater", "marine" , "marine.all",   
				#-> marine extant variables not significant
				#-> marine all near significant
				#->terrestrial not significant
				#->freshwater not significant
				
			#the following models ask: do turtles with/without neck retraction have a different mean labyrinth shapes?
			"incomplete_retr" , "full_retr",	#both non significant
			
			#the following models ask: do turtles with specific neck retractions have a different mean labyrinth shapes?
			"no_plane" , "vertical", "horizontal",  #all non significant
			
			"skull_box * skull_geometry + labyrinth_Csize + incomplete_retr",  #non-significant
			"skull_box * skull_geometry + labyrinth_Csize + full_retr", #non-significant
			"skull_box * skull_geometry + labyrinth_Csize + no_plane",  #non-significant
			"skull_box * skull_geometry + labyrinth_Csize + vertical", #non-significant
			"skull_box * skull_geometry + labyrinth_Csize + horizontal", #non-significant
				
			#check if ecological effects are redundant with skull size
			"skull_box + marine.all" , #marine.all becomes clearly non-significant, indicating the near-significant effect in bivarite models can be explained by skull size
			"labyrinth_Csize + marine.all" ,	#marine all remains near significnt significant
			"skull_geometry + marine.all" , #marine,all remains near significant
						
			#Further tests for effect of marine.all	
			"skull_box + skull_geometry + labyrinth_Csize + marine.all" , 
				#also insignificant, not further considered
			
			#are ecological effects important when included in the best model:
			"skull_box * skull_geometry + labyrinth_Csize + marine.all" ,
				# -> marine remains non significant
				
			#best model excludes ecological effects:
			"skull_box * skull_geometry + labyrinth_Csize"
				#-> all model parameters signifciant, brain case shape explains most of variance
			
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
					
				model.summaries[15]		
				
				#Print all coefficents to file
				capture.output(model.summaries, file = "Labyrinth_shape_model_summaries_incl_fossils.txt")
				
			##Altermative run with alternative tree
			
				gdf2 <- geomorph.data.frame( shape = GPA.labyrinth.all$coords[ ,, tree.temp$tip.label ] ,  
						phy = tree.temp.alternative ,
						marine = marine , freshwater = freshwater , terrestrial = terrestrial , marine.all = marine.all, 
						no_plane = no_plane, vertical = vertical, horizontal = horizontal,
						incomplete_retr = incomplete_retr, full_retr = full_retr,					 
						skull_length = log10( skull_length.temp ) , skull_width = log10(skull_width.temp) , skull_height = log10(skull_height.temp) , 
						skull_box = data.temp[ tree.temp$tip.label , "logV_mm3" ] ,  
						labyrinth_Csize = log10( labyrinth.Csize.all )[ tree.temp$tip.label ] ,
						skull_geometry =  skull_geometry.temp)
			
		
				procD.pgls.fits <- list()
					for( i in 1:length( models ) ) {
						procD.pgls.fits[[ i ]] <- procD.pgls( models[[ i ]] , phy = phy , SS.type = "II" , data = gdf2 )
					}
				
					##See summaries of procD.pgls results
						model.summaries <- lapply( procD.pgls.fits , summary )		
					
						#Print all coefficents to file
						capture.output(model.summaries, file = "Labyrinth_shape_model_summaries_incl_fossils_alternative_tree.txt")
					
		##Phylogenetic signal from Procrustes shape variables
		physig <- physignal(gdf$shape, phy = gdf$phy)
		summary(physig)
		plot(physig)
		
		