# Load packages

require( geomorph )
require( ape )
require( qpcR )
require( nlme )

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
			setwd( "INSERT DIRECTORY PATH"  )		
			tree <- read.nexus( "Dataset 5. cal3tree.calibrated.txt" )
			alternative.tree <- read.nexus("Dataset 6. mbltree.calibrated.txt")


		#Do GPA of labyrinth shape for extant taxa that are present in the phylogeny and have skull shape data
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
			
			#Tree for sensitivity analysis
			tree.temp.alternative <- drop.tip( alternative.tree , alternative.tree $tip.label[ ! alternative.tree $tip.label %in% names( labyrinth.Csize ) ] )

		#Prepare metadata that matches the labyrinth size data
			data.temp <- specimen.info[ specimen.info$Tree_names %in% tree.temp$tip.label , ]
				rownames( data.temp ) <-  data.temp$Tree_names
					data.temp <- data.temp[ tree.temp$tip.label , ]				
															 
		#Extract various things from metadata for analyses. 
		
			notwebbed <- data.temp$Forelimb == "not webbed - 0"
				names( notwebbed ) <- rownames( data.temp )
				
			minorlywebbed <- data.temp$Forelimb == "small - 1"
				names( minorlywebbed ) <- rownames( data.temp )
				
			moderatelywebbed <- data.temp$Forelimb == "extensive - 2"
				names( moderatelywebbed ) <- rownames( data.temp )
				
			extremelywebbed <- data.temp$Forelimb == "extensive - 3"
				names( extremelywebbed ) <- rownames( data.temp )

			flippered <- data.temp$Forelimb == "flipper - 4"
				names( flippered ) <- rownames( data.temp )
				
		##Habitat ecology	

			marine <-  data.temp$Habitat_general == "marine" 
				names( marine ) <- rownames( data.temp )	
					which(marine == TRUE)
				
			freshwater <-  data.temp$Habitat_general == "freshwater"  
				names( freshwater ) <- rownames( data.temp )
					which(freshwater == TRUE)
				
			terrestrial <-  data.temp$Habitat_general == "terrestrial"  
				names( terrestrial ) <- rownames( data.temp )
					which(terrestrial == TRUE)
				
		##Locomotor behaviour

			bottom_dwelling <-   data.temp$Habitual_habitat == "bottom_dwelling"  
				names( bottom_dwelling ) <- rownames( data.temp )	

			open_water <-  data.temp$Habitual_habitat == "open_water_swimming"  
				names( open_water ) <- rownames( data.temp )	

			digging <-  data.temp$Habitual_habitat == "digging" 
				names( digging ) <- rownames( data.temp )
				
			terrestrial_walking <-  data.temp$Habitual_habitat == "terrestrial_walking"   |  data.temp$Habitual_habitat == "digging" & !is.na( data.temp$Habitual_habitat ) 
				names( terrestrial_walking ) <- rownames( data.temp )
		
		##Neck retraction type categories
		
			vertical <-  data.temp $Retraction_type == "vertical" 
				names( vertical ) <- rownames( data.temp )
				
			horizontal <-  data.temp $Retraction_type == "sideways" 
				names( horizontal ) <- rownames( data.temp )

		##Neck retraction type categories
		
			full <-  data.temp $Retratction_ability == "full" 
				names( full ) <- rownames( data.temp )
				
			incomplete <-  data.temp $Retratction_ability == "incomplete" 
				names( incomplete ) <- rownames( data.temp )					
				
		##Size proxies
					
			skull_length.temp <- specimen.info$Skull_length_mm
				names(skull_length.temp) <- specimen.info$Tree_names
					skull_length.temp <- skull_length.temp[specimen.info$Tree_names %in% tree.temp$tip.label ]		
						skull_length.temp <- skull_length.temp[order(match(names(skull_length.temp), tree.temp$tip.label))  ]
														
			skull_width.temp <- specimen.info$Skull_width_mm
				names(skull_width.temp) <- specimen.info$Tree_names
					skull_width.temp <- skull_width.temp[specimen.info$Tree_names %in% tree.temp$tip.label ]
						skull_width.temp <- skull_width.temp[order(match(names(skull_width.temp), tree.temp$tip.label))  ]
					
			skull_height.temp <- specimen.info$Skull_height_mm
				names(skull_height.temp) <- specimen.info$Tree_names
					skull_height.temp <- skull_height.temp[specimen.info$Tree_names %in% tree.temp$tip.label ]
						skull_height.temp <- skull_height.temp[order(match(names(skull_height.temp), tree.temp$tip.label))  ]

		##Skull geometry proxy
		
				skull_geometry.temp <- skull_height.temp / skull_width.temp
				
					#check frequency distribution
					hist(skull_geometry.temp)
				
					#check if these make sense
					which(skull_geometry.temp[] == max(skull_geometry.temp))
					which(skull_geometry.temp[] == min(skull_geometry.temp))

		#Make a big data frame for analyses
		
			df <- data.frame( notwebbed, minorlywebbed, moderatelywebbed, extremelywebbed, flippered , 
							  marine , freshwater, terrestrial, 
							  bottom_dwelling , open_water , digging , terrestrial_walking, 
							  vertical, horizontal,
							  full, incomplete,
							  skull_length = log10( skull_length.temp ) , skull_width = log10(skull_width.temp), skull_height = log10(skull_height.temp), 
							  skull_box = data.temp[ tree.temp$tip.label , "logV_mm3" ] , 
							  labyrinth_Csize = log10( labyrinth.Csize )[ tree.temp$tip.label ] , 
							  skull_geometry =  skull_geometry.temp)	
							  				
				
				#Print table to file	
				write.table(df, file = "Labyrinth_CSize_explanatory_variables.csv", sep=" ")
				
		##In this script we're setting up all the combinations of explanatory variables for the right sizes of the models, and then running them all in a loop. This makes it easy to add regression models by extending the vector called "right.sides".
			right.sides <- c( "1" , 
			
			#following models ask can labyrinth size be predicted by skull size?
			"skull_length" , "skull_width", "skull_height", "skull_box" , 
			
			#can labyrinth size be predicted by the braincase aspect ratio?			
			"skull_geometry",
			
			#the following models ask: do terrestrial turtles have a different mean labyrinth size than non-terrestrial trutles?			
			"notwebbed", "minorlywebbed", "moderatelywebbed", "extremelywebbed", 
			"flippered" , "marine" , "freshwater", 
			"terrestrial", "bottom_dwelling" , "open_water" , "digging",  #"terrestrial_walking",
			
		#####the following models ask: do turtles with specific neck retraction types have a different mean labyrinth size ?
			"vertical", "horizontal",
					
		#####the following models ask: do turtles with the ability of neck retraction have a different mean labyrinth size ?
			"full", "incomplete",		
			
			#following models ask: is there an effect of braincase aspect ratio after skull size is taken into account? 
			"skull_box + skull_geometry", 
			
			#following models ask: do turtles of different webbing types have different labyrinth size means after skull size is accounted for?
			"skull_box + notwebbed",  "skull_box + minorlywebbed", "skull_box + moderatelywebbed", "skull_box + extremelywebbed", "skull_box + flippered",
			
			#following models ask: do turtles of different habitat preferences have different labyrinth size means after skull size is accounted for? 
			"skull_box + marine" , "skull_box + freshwater" ,"skull_box + terrestrial" ,
						
			#following models ask: do turtles of a specific locomotor behaviour have different labyrinth size means than other turtles after skull size is accounted for?
			#terrestrial_walking = terrestrial, thus muted 
			"skull_box + bottom_dwelling" , "skull_box + open_water" , "skull_box + digging" , #"skull_box + terrestrial_walking" ,
			
		#####the following models consider the neck retraction variables	
			"skull_box + vertical" , "skull_box + horizontal" ,
			"skull_box + full" , "skull_box + incomplete" ,			
			
			#following models ask: are there independent effect of skull size, braincase aspect ratio, and ecology
			"skull_box + skull_geometry + notwebbed",  "skull_box + skull_geometry + minorlywebbed", "skull_box + skull_geometry + moderatelywebbed", "skull_box + skull_geometry + extremelywebbed", "skull_box + skull_geometry + flippered",
			"skull_box + skull_geometry + marine" , "skull_box + skull_geometry + freshwater" ,"skull_box + skull_geometry + terrestrial" ,
			"skull_box + skull_geometry + bottom_dwelling" , "skull_box + skull_geometry + open_water" , "skull_box + skull_geometry + digging" , #"skull_box + skull_geometry + terrestrial_walking" ,	
			"skull_box + skull_geometry + vertical" , "skull_box + skull_geometry + horizontal" ,
			"skull_box + skull_geometry + full" , "skull_box + skull_geometry + incomplete" ,		 
						
			#the following models ask whether the included independent ecological effetcs account for any of the variance after the other effects were taken into account;
			#these specific effects were chosen because they were among best models in prior runs
			"skull_box + skull_geometry + freshwater + terrestrial",
			"skull_box + skull_geometry + open_water + notwebbed"
			
			)
			
				models <- paste( "labyrinth_Csize ~" , right.sides )
					models <- lapply( models , as.formula )

		#Run phylogenetic generalised least squares analyses (Grafen 1989) using Pagel's lambda correlation structure (Pagel 1999)
		
			pgls.fits <- list()
				for( i in 1:length( models ) ) {
					pgls.fits[[ i ]] <- gls( models[[ i ]] , correlation= corPagel( 0.9 , tree.temp ), data = df )
				}
				
		#Stores summary of models (used for later applications)
			model.summaries <- lapply( pgls.fits , summary )
							
		#Compute AICc values for all models. 
			AICc.values <- as.numeric(lapply( pgls.fits , AICc ))
				names(AICc.values) <- right.sides			
		#Store AICc, deltaAICc, relative likelihood and AICc weights into a dataframe, and order it by increasing AICc so that best model with lowest value is on top		
			AICc.delta.weigths <- data.frame(AICc.values, akaike.weights(AICc.values))
				AICc.delta.weigths.sorted <- AICc.delta.weigths[order(AICc.delta.weigths[,1]),]				
				#Check if by accident any of the models are duplicates, in which case it should print FALSE for the respective item
					length(AICc.delta.weigths.sorted$AICc.values) == length(unique(AICc.delta.weigths.sorted$AICc.values))
				
				##Add several columns of useful comparators to table
				#Calculate non-negligible AICc values based on 1/10th of AICc-weight of best model
					non.negligible <- rep("negligible", length(AICc.delta.weigths.sorted$weights))
						non.negligible.temp <- c(which(AICc.delta.weigths.sorted$weights > AICc.delta.weigths.sorted$weights[1]/10))
						non.negligible[non.negligible.temp] <- "non-negligible"
						non.negligible[1] <- "best model"				
					AICc.delta.weigths.sorted[,"AICc_weight_importance"] <- non.negligible
				#Add cumulative AICc weights
					AICc.delta.weigths.sorted[,"cum_AICc_weights"] <- cumsum(AICc.delta.weigths.sorted$weights)		
				#Add evidence ratios against best model; output value means 'best model is VALUE times likely to be the best model in AIC terms than is next best model'
						evidence_ratio <- c()
							for (i in 1:length(AICc.delta.weigths.sorted$weights)) {
								evidence_ratio[i] <- AICc.delta.weigths.sorted$weights[1]/AICc.delta.weigths.sorted$weights[i]
								}
					AICc.delta.weigths.sorted[,"evidence_ratio_bestmodel"] <- evidence_ratio
				#Add normalized preference probability against best model; ouput value means 'probability that best model is to be preferred over next best model is VALUE'
						norm_preference_probability	<- c()
							for (i in 1:length(AICc.delta.weigths.sorted$weights)) {
								norm_preference_probability[i] <- AICc.delta.weigths.sorted$weights[1]/(AICc.delta.weigths.sorted$weights[1] + AICc.delta.weigths.sorted$weights[i])
							}
					AICc.delta.weigths.sorted[,"norm_preference_prob_bestmodel"] <- norm_preference_probability
				#Add lambda value
						lambda.values <- c()
							for (i in 1:length(pgls.fits)) {
								lambda.values[i] <- as.numeric(pgls.fits[[i]]$modelStruct)
							}
						lambda.values.sorted <- lambda.values[order(AICc.delta.weigths[,1])]
					AICc.delta.weigths.sorted[,"lambda"] <- lambda.values.sorted									
				#Add R2 Nagelkerke 				
					#Define Log Likelihood of null model
					LogL_null <- pgls.fits[[1]]$logLik		
					#define order by which original model sequence has to be ordered to represent best, next-best, etc.
					besttoworst <- order(AICc.delta.weigths[,1])														
					#Extract all log liklihoods of all models (doesn't work with 'lapply( pgls.fits , logLik )' down the line because that also stored degrees of freedom etc, )
					LogL.all.models <- list()
					for( i in 1:length(pgls.fits) ) {
						LogL.all.models[i] <- pgls.fits[[i]]$logLik
						}					
					#Sort Log Likelihood according to AICc
					Log.all.models.sorted <- list()
					for( i in 1:length(LogL.all.models)) {
						Log.all.models.sorted[i] <- LogL.all.models[besttoworst[i]]
						}										
					#Calculate R2Nagerlke for all models
					Log.all.models.sorted.temp <- as.numeric(Log.all.models.sorted)
					R2.Nagelkerke.all.models <- list()
					for( i in 1:length(Log.all.models.sorted) ) {						
						R2.Nagelkerke.all.models[i] <- 1-exp(-(2/length(df$labyrinth_Csize)) * (Log.all.models.sorted.temp[i] - LogL_null) )						
					}					
					#Add R2 to table
					AICc.delta.weigths.sorted[,"R2"] <- as.numeric(R2.Nagelkerke.all.models)
					
				#Examine table
				AICc.delta.weigths.sorted
				
				#Print table to file	
				write.table(AICc.delta.weigths.sorted, file = "Labyrinth_CSize_Models_AICc_table.csv", sep=" ")			
									
				#Subset AICc table to only export best models for table in paper
				AICc.table.best.models <- AICc.delta.weigths.sorted[non.negligible.temp,]
					write.table(AICc.table.best.models, file = "Labyrinth_CSize_Models_AICc_table_best_models.csv", sep=" ")

		#Model coefficients are important extra information, but cannot be easily stored in same table, as each model has more than one coefficient
			#Extract model coefficients
				model.coefficients <- list()
					for (i in 1:length(models) ) {
						model.coefficients[[i]] <- coef(model.summaries[[i]])				
					}						
			#Order model coefficients according to increasing AICc scores of models	
				besttoworst <- order(AICc.delta.weigths[,1]) #order by which original model sequence has to be ordered to represent best, next-best, etc.
					
				model.coefficients.ordered <- list()
				for( i in 1:length(model.coefficients)) {
					model.coefficients.ordered[i] <- model.coefficients[besttoworst[i]]
					}
			#Print all coefficents to file
				capture.output(model.coefficients.ordered, file = "Labyrinth_CSize_Coefficients_all_models_ordered.txt")
			#Extract coefficients for non-negigible models, as these will have to be reported and it's easier to oversee than all models' coefficients		
				coefficients.non.negigible.models <- model.coefficients.ordered[non.negligible.temp]
					capture.output(coefficients.non.negigible.models, file = "Labyrinth_CSize_Coefficients_non_negigible_models.txt")
				
				
	####Sensitivity analysis with different tree		
	
		#PGLS analysis		
			pgls.fits2 <- list()
			for( i in 1:length( models ) ) {
				pgls.fits2[[ i ]] <- gls( models[[ i ]] , correlation= corPagel( 0.9 , tree.temp.alternative ), data = df )
			}			
			
			#Stores summary of models (used for later applications)
			model.summaries2 <- lapply( pgls.fits2 , summary )		

			#Compute AICc values for all models. 
			AICc.values <- as.numeric(lapply( pgls.fits2 , AICc ))
				names(AICc.values) <- right.sides			
		
		#Store AICc, deltaAICc, relative likelihood and AICc weights into a dataframe, and order it by increasing AICc so that best model with lowest value is on top		
			AICc.delta.weigths <- data.frame(AICc.values, akaike.weights(AICc.values))
				AICc.delta.weigths.sorted <- AICc.delta.weigths[order(AICc.delta.weigths[,1]),]				
				#Check if by accident any of the models are duplicates, in which case it should print FALSE for the respective item
					length(AICc.delta.weigths.sorted$AICc.values) == length(unique(AICc.delta.weigths.sorted$AICc.values))
				
				##Add several columns of useful comparators to table
				#Calculate non-negligible AICc values based on 1/10th of AICc-weight of best model
					non.negligible <- rep("negligible", length(AICc.delta.weigths.sorted$weights))
						non.negligible.temp <- c(which(AICc.delta.weigths.sorted$weights > AICc.delta.weigths.sorted$weights[1]/10))
						non.negligible[non.negligible.temp] <- "non-negligible"
						non.negligible[1] <- "best model"				
					AICc.delta.weigths.sorted[,"AICc_weight_importance"] <- non.negligible
				#Add cumulative AICc weights
					AICc.delta.weigths.sorted[,"cum_AICc_weights"] <- cumsum(AICc.delta.weigths.sorted$weights)		
				#Add evidence ratios against best model; output value means 'best model is VALUE times likely to be the best model in AIC terms than is next best model'
						evidence_ratio <- c()
							for (i in 1:length(AICc.delta.weigths.sorted$weights)) {
								evidence_ratio[i] <- AICc.delta.weigths.sorted$weights[1]/AICc.delta.weigths.sorted$weights[i]
								}
					AICc.delta.weigths.sorted[,"evidence_ratio_bestmodel"] <- evidence_ratio
				#Add normalized preference probability against best model; ouput value means 'probability that best model is to be preferred over next best model is VALUE'
						norm_preference_probability	<- c()
							for (i in 1:length(AICc.delta.weigths.sorted$weights)) {
								norm_preference_probability[i] <- AICc.delta.weigths.sorted$weights[1]/(AICc.delta.weigths.sorted$weights[1] + AICc.delta.weigths.sorted$weights[i])
							}
					AICc.delta.weigths.sorted[,"norm_preference_prob_bestmodel"] <- norm_preference_probability
				#Add lambda value
						lambda.values <- c()
							for (i in 1:length(pgls.fits)) {
								lambda.values[i] <- as.numeric(pgls.fits[[i]]$modelStruct)
							}
						lambda.values.sorted <- lambda.values[order(AICc.delta.weigths[,1])]
					AICc.delta.weigths.sorted[,"lambda"] <- lambda.values.sorted									
				#Add R2 Nagelkerke 				
					#Define Log Likelihood of null model
					LogL_null <- pgls.fits2[[1]]$logLik		
					#define order by which original model sequence has to be ordered to represent best, next-best, etc.
					besttoworst <- order(AICc.delta.weigths[,1])														
					#Extract all log liklihoods of all models (doesn't work with 'lapply( pgls.fits , logLik )' down the line because that also stored degrees of freedom etc, )
					LogL.all.models <- list()
					for( i in 1:length(pgls.fits2) ) {
						LogL.all.models[i] <- pgls.fits2[[i]]$logLik
						}					
					#Sort Log Likelihood according to AICc
					Log.all.models.sorted <- list()
					for( i in 1:length(LogL.all.models)) {
						Log.all.models.sorted[i] <- LogL.all.models[besttoworst[i]]
						}										
					#Calculate R2Nagerlke for all models
					Log.all.models.sorted.temp <- as.numeric(Log.all.models.sorted)
					R2.Nagelkerke.all.models <- list()
					for( i in 1:length(Log.all.models.sorted) ) {						
						R2.Nagelkerke.all.models[i] <- 1-exp(-(2/length(df$labyrinth_Csize)) * (Log.all.models.sorted.temp[i] - LogL_null) )						
					}					
					#Add R2 to table
					AICc.delta.weigths.sorted[,"R2"] <- as.numeric(R2.Nagelkerke.all.models)
					
				#Examine table
				AICc.delta.weigths.sorted
				
				#Print table to file	
				write.table(AICc.delta.weigths.sorted, file = "Labyrinth_CSize_Models_AICc_table_alternative_tree.csv", sep=" ")			
									
				#Subset AICc table to only export best models for table in paper
				AICc.table.best.models <- AICc.delta.weigths.sorted[non.negligible.temp,]
					write.table(AICc.table.best.models, file = "Labyrinth_CSize_Models_AICc_table_best_models_alternative_tree.csv", sep=" ")

		#Model coefficients are important extra information, but cannot be easily stored in same table, as each model has more than one coefficient
			#Extract model coefficients
				model.coefficients <- list()
					for (i in 1:length(models) ) {
						model.coefficients[[i]] <- coef(model.summaries2[[i]])				
					}						
			#Order model coefficients according to increasing AICc scores of models	
				besttoworst <- order(AICc.delta.weigths[,1]) #order by which original model sequence has to be ordered to represent best, next-best, etc.
					
				model.coefficients.ordered <- list()
				for( i in 1:length(model.coefficients)) {
					model.coefficients.ordered[i] <- model.coefficients[besttoworst[i]]
					}
			#Print all coefficents to file
				capture.output(model.coefficients.ordered, file = "Labyrinth_CSize_Coefficients_all_models_ordered_alternative_tree.txt")
			#Extract coefficients for non-negigible models, as these will have to be reported and it's easier to oversee than all models' coefficients		
				coefficients.non.negigible.models <- model.coefficients.ordered[non.negligible.temp]
					capture.output(coefficients.non.negigible.models, file = "Labyrinth_CSize_Coefficients_non_negigible_models_alternative_tree.txt")
				