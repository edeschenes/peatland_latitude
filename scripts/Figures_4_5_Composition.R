
#Load packages
require(eulerr)
require(vegan) 
require(FD)
require(svglite)

set_theme(base = theme_classic()) 

Sys.setenv(LANGUAGE="en")

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

#Select only environmental variables used in the composition analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center). The unstandardized data will be used to obtain predicted values from the models to be used in the plots (Figure 1). 
env[,c(2:7)] <- scale(env[,c(2:7)], center = TRUE, scale = TRUE)
env$Habitat <- factor(env$Habitat)

##############################################################################
####RDA of taxonomic composition as a function of environmental variables ####
##############################################################################

#Hellinger transformation of the site x species abundance matrix
vas.hel <- decostand(vascu, method = 'hellinger')
bry.hel <- decostand(bryo, method = 'hellinger')

# RDA of hellinger transformed vascular taxonomic composition as a function of standardized environmental variables
(vas.rda <- rda(vas.hel ~ ., env))
summary(vas.rda)	# Scaling 2 (default) 
anova(vas.rda) # is the model significant?

# RDA of hellinger transformed bryophyte taxonomic composition as a function of standardized environmental variables
(bry.rda <- rda(bry.hel ~ ., env)) 
summary(bry.rda)	# Scaling 2 (default) 
anova(bry.rda) # is the model significant?

#The RDA models just performed will be used for the variation partitioning. 

#Define the groups for environmental variables: bioclimatic, spatial, local
bioclim <- env[, c(4:5)]
space <- env[, c(2:3)]
local <- env[, c(1,6,7)]


##################################################################
#### Variation partitioning of vascular taxonomic composition ####
##################################################################

(spe.part.vas <- varpart(vascu, bioclim,space, local, transfo= 'hel'))

#Significance of each testable fraction of the variation partitioning (unique fraction)
(frac.1 <- anova(rda(vas.hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(vas.hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(vas.hel, local, cbind.data.frame(bioclim, space)), step=10000))

#Plot of variation partitioning Venn diagram
plot(spe.part.vas, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim', 'Space', 'Local'))

#Assign to a Euler object the proportion of variation explained by each group of variables and their shared proportions
venn_vasc_tax <- euler(c(A = 1,          
                     B = 1.4,
                     C = 12.4,
                     "A&B" = 2.1,
                     'A&C'= 0.5,
                     'B&C' = 0.5,
                     'A&B&C' = 0))

#Plot the Venn diagram with circle sizes proportional to the proportion of variation explained
p_venn_vasc_tax <- plot(venn_vasc_tax, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)

##############################################################
#### Variation partitioning of moss taxonomic composition ####
##############################################################

(spe.part.bry <- varpart(bryo, bioclim,space, local,transfo= 'hel'))

#Significance of each testable fraction of the variation partitioning (unique fraction)
(frac.1 <- anova(rda(bry.hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(bry.hel, local, cbind.data.frame(bioclim, space)), step=10000))
(frac.3 <- anova(rda(bry.hel, space, cbind.data.frame(bioclim, local)), step=10000))

#Plot of variation partitioning Venn diagram
plot(spe.part.bry, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim', 'Space', 'Local'))

#Assign to a Euler object the proportion of variation explained by each group of variables and their shared proportions
venn_bryo_tax <- euler(c(A = 0.6,          # Draw pairwise venn diagram
                         B = 0.9,
                         C = 6.9,
                         "A&B" = 1.9,
                         'A&C'= 0.3,
                         'B&C' = 0.3,
                         'A&B&C' = 0.5))

#Plot the Venn diagram with circle sizes proportional to the proportion of variation explained
p_venn_bryo_tax <- plot(venn_bryo_tax, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)

###############################################################################
####RDA of functional composition as a function of environmental variables ####
###############################################################################
#RDA is conducted on the site x trait class frequencies matrix.
#Step 1. Calculate trait class frequencies with function functcomp (CWM)
#Step 2. PCA on Hellinger transformed site x trait class frequencies matrix

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

#Select only traits used for calculation of FDis
traits_vasc <- dplyr::select(traits_vasc, Port, Height, Seed, SLA)

#Order site x vascular abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species. 
vascu <-vascu[,order(colnames(vascu))]
#Create site x vascular abundance matrix
vascu<- as.matrix(vascu)

#### Calculate vascular CWM ####
cwm_vascu <- functcomp(traits_vasc, vascu, CWM.type = 'all')# returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_vascu) <- c('Herbaceous', 'Shrub', 'Tree','Height 1', 'Height 2', 'Height 3', 'Seed 1', 'Seed 2', 'Seed 3', 'SLA 1', 'SLA 2', 'SLA 3')

####################################
#### Calculate moss species CWM ####
####################################

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

#Select only traits used for calculation of FDis
traits_bryo <- dplyr::select(traits_bryo, Bryophyte_group, Life_form, Shoot_length, Sexual_condition, Seta_length, Peristome_type, Spore_size, Tomentum)
sapply(traits_bryo, class)

#Order site x moss abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species.
bryo <-bryo[,order(colnames(bryo))]

#Create site x vascular abundance matrix
bryo<- as.matrix(bryo)

#### Calculate CWM ####
cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_bryo) <- c('Acrocarpous','Pleurocarpous','Sphagnum','Tuft','Turf','Weft','Shoot 1','Shoot 2','Shoot 3','Dioicous','Monoicous','Seta 0','Seta 2','Seta 3','Seta 4','Diplolepidous alternate','Diplolepidous perfect','Haplolepideous','Nematodontus','No peristome','Spore 1','Spore 2','Spore 3','No_Tomentum', 'Tomentum')

cwm_bryo <- subset(cwm_bryo, select = -c(No_Tomentum, Monoicous))

##############################################################################
####RDA of taxonomic composition as a function of environmental variables ####
##############################################################################

#Hellinger transformation of the site x trait class frequencies matrix
traits_vasc_hel <- sqrt(cwm_vascu)
traits_bryo_hel <- sqrt(cwm_bryo)


#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center).
env.scale <- env
env.scale[,c(2:7)] <- scale(env.scale[,c(2:7)], center = TRUE, scale = TRUE)

# RDA of hellinger transformed vascular CWM as a function of standardized environmental variables
(vas.funct.rda <- rda(traits_vasc_hel ~ ., env.scale))
summary(vas.funct.rda)	# Scaling 2 (default) 
anova(vas.funct.rda) # is the model significant?

# RDA of hellinger transformed bryophyte CWM as a function of standardized environmental variables
(bry.funct.rda <- rda(traits_bryo_hel ~ ., env.scale)) 
summary(bry.funct.rda)	# Scaling 2 (default) 
anova(bry.funct.rda) # is the model significant?

#The RDA models just performed will be used for the variation partitioning. 

#Define the groups for environmental variables: bioclimatic, spatial, local
bioclim <- env[, c(4:5)]
space <- env[, c(2:3)]
local <- env[, c(1,6,7)]


###################################################################
#### Variation partitioning of vascular functional composition ####
###################################################################

(spe.part.vas <- varpart(traits_vasc_hel, bioclim,space, local))

#Significance of each testable fraction of the variation partitioning (unique fraction)
(frac.1 <- anova(rda(traits_vasc_hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(traits_vasc_hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(traits_vasc_hel, local, cbind.data.frame(bioclim, space)), step=10000))

#Plot of variation partitioning Venn diagram
plot(spe.part.vas, digits = 2, bg = c("red", "blue", 'green','yellow'),Xnames = c('Bioclimatique', 'Space', 'Local'))

#Assign to a Euler object the proportion of variation explained by each group of variables and their shared proportions
venn_vasc_funct <- euler(c(A = 0.5,          # Draw pairwise venn diagram
                           B = 1.5,
                           C = 13.5,
                           "A&B" = 2.8,
                           'A&C'= 0.1,
                           'B&C' = 0.5,
                           'A&B&C' =0))

#Plot the Venn diagram with circle sizes proportional to the proportion of variation explained
p_venn_vasc_funct<- plot(venn_vasc_funct, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)


###############################################################
#### Variation partitioning of moss functional composition ####
###############################################################

(spe.part.bry <- varpart(traits_bryo_hel, bioclim,space, local))

#Significance of each testable fraction of the variation partitioning (unique fraction)
(frac.1 <- anova(rda(traits_bryo_hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(traits_bryo_hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(traits_bryo_hel, local, cbind.data.frame(bioclim, space)), step=10000))

#Plot of variation partitioning Venn diagram
plot(spe.part.bry, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim','Space', 'Local'))

#Assign to a Euler object the proportion of variation explained by each group of variables and their shared proportions
venn_bryo_funct <- euler(c(A = 2.1,          
                           B = 1.6,
                           C = 6.5,
                           "A&B" = 1.8,
                           'A&C'= 1.1,
                           'B&C' = 0.7,
                           'A&B&C' = 0))

#Plot the Venn diagram with circle sizes proportional to the proportion of variation explained
p_venn_bryo_funct <- plot(venn_bryo_funct, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)

#######################################################
####Download SVG file for modification in Inkscape ####
#######################################################
svglite("figure4_varpart.svg")
gridExtra::grid.arrange(p_venn_vasc_tax, p_venn_bryo_tax,p_venn_vasc_funct, p_venn_bryo_funct, ncol=2)
dev.off()

##################################
#### RDA OF FUNCTIONAL TRAITS ####
##################################

#Assign colnames for plotting
colnames(env) <- c("Habitat", "Latitude",'Longitude', 'Annual temperature', 'Annual precipitations', 'Thickness', 'Surface water')

#Assign habitat as factor for plotting
env$Habitat=as.factor(env$Habitat)

#Prep for RDA plot
par(mfrow=c(1,1))
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)

#################################################################################
#### Plot of vascular species functional traits as a function of enviro variables
#################################################################################
#Both RDA plots will be combined and cleaned in Inkscape. 

summary(vas.funct.rda)

svglite("figure5_RDA_vasc.svg")
plot(vas.funct.rda, scaling=2, main="Triplot RDA - vascular traits", type="none", xlab=c("RDA1 (88.45%)"), ylab=c("RDA2 (4.74%)"), xlim=c(-1, 1), ylim=c(-1.4,1.4), cex =1)
with(env, points(vas.funct.rda, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(vas.funct.rda, display="species", choices=c(1), scaling=2)*2,
       scores(vas.funct.rda, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(vas.funct.rda, display="species", choices=c(1), scaling=2)*2.4,
     scores(vas.funct.rda, display="species", choices=c(2), scaling=2)*2.4,
     labels=rownames(scores(vas.funct.rda, display="species", scaling=2)),
     col="black", cex=1)   
arrows(0,0,
       scores(vas.funct.rda, display="bp", choices=c(1), scaling=2)*2,
       scores(vas.funct.rda, display="bp", choices=c(2), scaling=2)*2,
       col="red",length = 0.1)
text(scores(vas.funct.rda, display="bp", choices=c(1), scaling=2)*2,
     scores(vas.funct.rda, display="bp", choices=c(2), scaling=2)*2,
     labels=rownames(scores(vas.funct.rda, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1) 
dev.off()


#################################################################################
#### Plot of bryophyte species functional traits as a function of enviro variables
#################################################################################
summary(bry.funct.rda)

svglite("figure5_RDA_bryo.svg")
plot(bry.funct.rda, scaling=2, main="Triplot RDA - bryo traits", type="none", xlab=c("RDA1 (55.01%)"), ylab=c("RDA2 (29.03%)"), cex = 1,ylim=c(-1.3,1.3))
with(env, points(bry.funct.rda, display = "sites", col = colvec[Habitat],
                 scaling = 2, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec,pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(bry.funct.rda, display="species", choices=c(1), scaling=2)*2,
       scores(bry.funct.rda, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(bry.funct.rda, display="species", choices=c(1), scaling=2)*2.1,
     scores(bry.funct.rda, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(bry.funct.rda, display="species", scaling=2)),
     col="black", cex=1)   
arrows(0,0,
       scores(bry.funct.rda, display="bp", choices=c(1), scaling=2)*1.5,
       scores(bry.funct.rda, display="bp", choices=c(2), scaling=2)*1.5,
       col="red",length = 0.1)
text(scores(bry.funct.rda, display="bp", choices=c(1), scaling=2)*1.5+0.05,
     scores(bry.funct.rda, display="bp", choices=c(2), scaling=2)*1.5+0.05,
     labels=rownames(scores(bry.funct.rda, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)
dev.off()
