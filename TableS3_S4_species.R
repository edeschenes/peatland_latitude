###Descriptive tables S3 and S4. 
##Load packages


#Load data
#Site x environmental variables matrix
sites <- read.csv("C:/Users/edeschen/Documents/MAITRISE/Scripts/Data_repository/sites_data.csv")
#Site x species abundance matrix for vascular species
vascu <- read.csv("C:/Users/edeschen/Documents/MAITRISE/Scripts/Data_repository/vascu_data.csv")
#Site x species abundance matrix for moss species
bryo <- read.csv("C:/Users/edeschen/Documents/MAITRISE/Scripts/Data_repository/bryo_data.csv")

#Assign plotID to row names and remove plotID as a column
vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)
sites <- data.frame(sites, row.names = sites$ID_Quadrat)

env <- dplyr::select(sites,Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

##
########################################
#### Calculate vascular species CWM ####
########################################

#Open vascular trait data 
traits_vasc <- read.csv("C:/Users/edeschen/Documents/MAITRISE/Scripts/Data_repository/traits_vasc_data.csv",sep = ";", header= TRUE)

#Assign columns as factor
traits_vasc <- as.data.frame(unclass(traits_vasc), stringsAsFactors = TRUE)

#Assign ordinal columns as ordered 
traits_vasc$SLA <- factor(traits_vasc$SLA, order = TRUE, 
                          levels = c("SLA 1", "SLA 2", "SLA 3"))
traits_vasc$Height <- factor(traits_vasc$Height, order = TRUE, 
                             levels = c("Height 1", "Height 2", "Height 3"))
traits_vasc$Seed <- factor(traits_vasc$Seed, order = TRUE, 
                           levels = c("Seed 1", "Seed 2", "Seed 3"))

#Assign species code as row names
traits_vasc <- data.frame(traits_vasc, row.names = traits_vasc$Code)

#Order by row names
traits_vasc <- traits_vasc[order(traits_vasc$Code),]

#Order site x vascular abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species. 
vascu <-vascu[,order(colnames(vascu))]
#Create site x vascular abundance matrix
vascu<- as.matrix(vascu)


#Open moss trait data 
traits_bryo <- read.csv("C:/Users/edeschen/Documents/MAITRISE/Scripts/Data_repository/traits_moss_data.csv", sep = ";", header= TRUE)

#Assign columns as factor
traits_bryo <- as.data.frame(unclass(traits_bryo), stringsAsFactors = TRUE)

#Assign ordinal columns as ordered 
traits_bryo$Shoot_length <- factor(traits_bryo$Shoot_length, order = TRUE, 
                                   levels = c("Shoot 1", "Shoot 2", "Shoot 3"))
traits_bryo$Seta_length <- factor(traits_bryo$Seta_length, order = TRUE, 
                                  levels = c('Seta 0', "Seta 2", "Seta 3", 'Seta 4'))
traits_bryo$Spore_size <- factor(traits_bryo$Spore_size, order = TRUE, 
                                 levels = c("Spore 1", "Spore 2", "Spore 3"))

#Assign species code as row names
traits_bryo <- data.frame(traits_bryo, row.names = traits_bryo$Code)

#Order by row names
traits_bryo <- traits_bryo[order(traits_bryo$Code),]

#Order site x moss abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species.
bryo <-bryo[,order(colnames(bryo))]

#Create site x vascular abundance matrix
bryo<- as.matrix(bryo)

###############################################################################

##Select bogs and subset all dataframes
bog <- sites[sites[, 'Habitat'] == 'Bog', ] 
quad_bog <- unique(bog$ID_Quadrat)

#Subset site x species matrices
bog_bry <- bryo[row.names(bryo) %in% quad_bog,]
bog_vas <- vascu[row.names(vascu) %in% quad_bog,]

#Remove species that are not present in bogs
bog_vas <- bog_vas[,colSums(bog_vas != 0) > 0]
bog_bry <- bog_bry[,colSums(bog_bry != 0) > 0]

#Subset site x environmental variables matrix
env_bog<- env[row.names(env) %in% quad_bog,] #169 sites
env_bog <- subset(env_bog, select = c(-Habitat))
env_bog <- as.data.frame(scale(env_bog, center = TRUE, scale = TRUE))

##Select fens and subset all dataframes
fen <- sites[sites[, 'Habitat'] == 'Fen', ] 
quad_fen <- unique(fen$ID_Quadrat)

#Subset site x species matrices
fen_bry <- bryo[row.names(bryo) %in% quad_fen,]
fen_vas <- vascu[row.names(vascu) %in% quad_fen,]

#Subset site x environmental variables matrix
env_fen<- env[row.names(env) %in% quad_fen,] #231 sites
env_fen <- subset(env_fen, select = c(-Habitat))
env_fen <- as.data.frame(scale(env_fen, center = TRUE, scale = TRUE))


#Calculate frequencies of species within bogs
spc.pres.bog.vasc <- as.data.frame(apply(bog_vas>0,2,sum))
spc.pres.bog.bry <- as.data.frame(apply(bog_bry>0,2,sum))

#Merge species frequencies with functional trait data
traits_vasc <- merge(traits_vasc,spc.pres.bog.vasc,by="row.names",all.x=TRUE)
traits_bryo <- merge(traits_bryo,spc.pres.bog.bry,by="row.names",all.x=TRUE)

#Remove duplicate row names and assign row names
traits_vasc <- data.frame(traits_vasc[,-1], row.names = traits_vasc$Row.names)
traits_bryo <- data.frame(traits_bryo[,-1], row.names = traits_bryo$Row.names)

#Calculate frequencies of species within fens
spc.pres.fen.vasc <- as.data.frame(apply(fen_vas>0,2,sum))
spc.pres.fen.bry <- as.data.frame(apply(fen_bry>0,2,sum))

#Merge species frequencies with functional trait data
traits_vasc <- merge(traits_vasc,spc.pres.fen.vasc,by="row.names",all.x=TRUE)
traits_bryo <- merge(traits_bryo,spc.pres.fen.bry,by="row.names",all.x=TRUE)

#Remove duplicate row names and assign row names
traits_vasc <- data.frame(traits_vasc[,-1], row.names = traits_vasc$Row.names)
traits_bryo <- data.frame(traits_bryo[,-1], row.names = traits_bryo$Row.names)

##Change colnames
colnames(traits_vasc) <- c('Species','Code','Life_form','Max_height_cm','Seed_weight_mg','SLA_num','Seed','Height','SLA','Frenquency_bog','Frequency_fen')
colnames(traits_bryo) <- c('Species','Code','Bryophyte_group','Life_form','Shoot_length_cm','Shoot_length','Sexual_condition','Seta_length','Peristome_type','Spore_size_um','Spore_size','Tomentum', 'Frequency_bog','Frequency_fen')

#Export tables
write.csv(traits_vasc, 'table_S3_vasc.csv')
write.csv(traits_bryo, 'table_S4_bryo.csv')

