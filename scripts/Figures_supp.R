#Scripts to create figures S2, S3 and S4. 
library(ggcorrplot)
require(GGally)
require(vegan) 
require(FD)
require(svglite)
require(gridExtra)
library(reshape2)
require(ggpmisc)

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

#Select only environmental variables used in the supplementary analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

#Assign habitat as factor
env$Habitat=as.factor(env$Habitat)

#################################################
#### Seperate bogs and fens in all dataframes####
#################################################

####Bogs ####
bog <- sites[sites[, 'Habitat'] == 'Bog', ] #select bogs 
quad_bog <- unique(bog$ID_Quadrat) #assign bog quadrat IDs to object

bog_vas <- vascu[row.names(vascu) %in% quad_bog,]#select bogs in site x vascular abundance matrix
bog_bry <- bryo[row.names(bryo) %in% quad_bog,] #select bogs in site x moss abundance matrix

#Remove species that are not present in bogs
bog_vas <- bog_vas[,colSums(bog_vas != 0) > 0]
bog_bry <- bog_bry[,colSums(bog_bry != 0) > 0]

#Create environmental dataframe with only bogs
env_bog<- env[row.names(env) %in% quad_bog,] #select bogs in site x enviro matrix
#Scale and center environmental variables
env_bog[,c(2:7)] <- as.data.frame(scale(env_bog[,c(2:7)], center = TRUE, scale = TRUE))

####Fens ####
fen <- sites[sites[, 'Habitat'] == 'Fen', ] #select fens
quad_fen <- unique(fen$ID_Quadrat) #assign fen quadrat IDs to object

fen_vas <- vascu[row.names(vascu) %in% quad_fen,] #select fens in site x vascular abundance matrix
fen_bry <- bryo[row.names(bryo) %in% quad_fen,] #select fens in site x moss abundance matrix

#Create environmental dataframe with only bogs
env_fen<- env[row.names(env) %in% quad_fen,]#select fens in site x enviro matrix
#Scale and center environmental variables
env_fen[,c(2:7)] <- as.data.frame(scale(env_fen[,c(2:7)], center = TRUE, scale = TRUE))


##################################################
#### Figure S1 - PCA of bioclimatic variables ####
##################################################

#Select all bioclimatic variables
env1 <- dplyr::select(sites,bioclim_1, bioclim_2, bioclim_3, bioclim_4, bioclim_5, bioclim_6, bioclim_7, bioclim_8, bioclim_9, bioclim_10, bioclim_11, bioclim_12, bioclim_13, bioclim_14, bioclim_15, bioclim_16, bioclim_17, bioclim_18, bioclim_19)

#Run PCA of bioclimatic variables
env.pca <- rda(env1, scale = TRUE)
summary(env.pca) # Default scaling 2

#Prep for PCA plot
par(mfrow=c(1,1))
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)
sites$Habitat=as.factor(sites$Habitat)

#Plot of PCA of bioclimatic variables
svglite("figure_S1_pca_bioclim.svg")
plot(env.pca, scaling=2, type="none", xlab=c("PCA1 (47.73%)"), ylab=c("PCA2 (31.03%)"),xlim=c(-2,3), ylim=c(-1.5, 2))
points(scores(env.pca, display="sites", choices=c(1,2), scaling=2),pch=21, col="black", bg="black", cex=0.7)
with(sites, points(env.pca, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(sites, legend("topleft", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(env.pca, display="species", choices=c(1), scaling=2),
       scores(env.pca, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(env.pca, display="species", choices=c(1), scaling=2),
     scores(env.pca, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(env.pca, display="species", scaling=2)),
     col="black", cex=1)   
dev.off()

#################################################################
#### Figure S2 - Correlation plot of environmental variables ####
#################################################################

#Select environmental variables in bogs
env_bog<- dplyr::select(env_bog, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

#Change column names for plotting
colnames(env_bog) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water')

#Select environmental variables in fens
env_fen<- dplyr::select(env_fen, Latitude,Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

#Change column names for plotting
colnames(env_fen) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water')

#Calculate Pearson correlations between each pair of environmental variables
r_bog <- cor(env_bog, use="complete.obs")
r_fen <- cor(env_fen, use="complete.obs")

#Plot pairwise plot of Pearson correlations between environmental variables 
cor_bog <- ggcorr(env_bog, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)
cor_fen <- ggcorr(env_fen, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)

#######################################################
####Download SVG file for modification in Inkscape ####
#######################################################
svglite("figure_S2_corplot.svg")
gridExtra::grid.arrange(cor_bog, cor_fen, ncol=2)
dev.off()

##############################################################
#### Figure S3 - Environmental variables in bogs and fens ####
##############################################################
#Select environmental variables
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

####Plot each variable as a functional of latitude, separately for bogs and fens 

p1<-ggplot(data= env, aes(Habitat, Latitude))+
  stat_boxplot( aes(Habitat, Latitude), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Latitude),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Latitude (degrees)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Latitude~Habitat, data=env) #run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

p2 <-ggplot(data= env, aes(Habitat, Longitude))+
  stat_boxplot( aes(Habitat, Longitude), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Longitude),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Longitude (degrees)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Longitude~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

p3 <-ggplot(data= env, aes(Habitat, bioclim_1))+
  stat_boxplot( aes(Habitat, bioclim_1), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, bioclim_1),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Annual temperature (?C)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(bioclim_1~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

p4<-ggplot(data= env, aes(Habitat, bioclim_12))+
  stat_boxplot( aes(Habitat, bioclim_12), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, bioclim_12),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Annual precipitation (mm)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(bioclim_12~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

p5<-ggplot(data= env, aes(Habitat, Thickness))+
  stat_boxplot( aes(Habitat, Thickness), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Thickness),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Peat thickness (cm)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Thickness~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

p6<-ggplot(data= env, aes(Habitat, Surface_water))+
  stat_boxplot( aes(Habitat, Surface_water), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Surface_water),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3)+
  labs(y="Surface water (%)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Surface_water~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

#######################################################
####Download SVG file for modification in Inkscape ####
#######################################################
svglite("figure_S3_env.svg", width = 14, height = 7)
grid.arrange(p1,p3,p5,p2,p4,p6,ncol=3)
dev.off()

############################################
####Figure S4. PCA of functional traits ####
############################################
#PCA is conducted on the site x trait class frequencies matrix.
#Step 1. Calculate trait class frequencies with function functcomp
#Step 2. PCA on Hellinger transformed site x trait class frequencies matrix

####

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

#Select only traits used for calculation of CWM
traits_bryo <- dplyr::select(traits_bryo, Bryophyte_group, Life_form, Shoot_length, Sexual_condition, Seta_length, Peristome_type, Spore_size, Tomentum)

#Order site x moss abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species.
bryo <-bryo[,order(colnames(bryo))]

#Create site x vascular abundance matrix
bryo<- as.matrix(bryo)

#### Calculate CWM ####
cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_bryo) <- c('Acrocarpous','Pleurocarpous','Sphagnum','Tuft','Turf','Weft','Shoot 1','Shoot 2','Shoot 3','Dioicous','Monoicous','Seta 0','Seta 2','Seta 3','Seta 4','Diplolepidous alternate','Diplolepidous perfect','Haplolepideous','Nematodontus','No peristome','Spore 1','Spore 2','Spore 3','No_Tomentum', 'Tomentum')

cwm_bryo <- subset(cwm_bryo, select = -c(No_Tomentum, Monoicous))


#Prep for PCA plot
par(mfrow=c(1,1))
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)
sites$Habitat=as.factor(sites$Habitat)


####################################################
#### Plot of vascular species functional traits ####
####################################################
traits_vasc_hel <- sqrt(cwm_vascu) #Hellinger transformation on site x trait class frequencies
vas.funct.pca <- rda(traits_vasc_hel) #perform PCA for vascular traits
summary(vas.funct.pca)

#Plot
svglite("figure_S4_pca_vasc.svg")
plot(vas.funct.pca, scaling=2, main="PCA - vascular traits", type="none", xlab=c("PCA1 (56.02%)"), ylab=c("PCA2 (12.92%)"), xlim=c(-1.3, 1.3), ylim=c(-1,1), cex =1)
with(env, points(vas.funct.pca, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(vas.funct.pca, display="species", choices=c(1), scaling=2),
       scores(vas.funct.pca, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(vas.funct.pca, display="species", choices=c(1), scaling=2)*1.1,
     scores(vas.funct.pca, display="species", choices=c(2), scaling=2)*1.1,
     labels=rownames(scores(vas.funct.pca, display="species", scaling=2)),
     col="black", cex=1)   
dev.off()


################################################
#### Plot of moss species functional traits ####
################################################
traits_bryo_hel <- sqrt(cwm_bryo)#Hellinger transformation on site x trait class frequencies
bry.funct.pca <- rda(cwm_bryo) #perform PCA for moss traits
summary(bry.funct.pca)

#Plot
svglite("figure_S4_pca_bryo.svg")
plot(bry.funct.pca, scaling=2, main="PCA - bryo traits", type="none", xlab=c("PCA1 (42.45%)"), ylab=c("PCA2 (28.58%)"), xlim=c(-1.3, 1.3), ylim=c(-1.6,1.6), cex = 1)
with(env, points(bry.funct.pca, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec,pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(bry.funct.pca, display="species", choices=c(1), scaling=2),
       scores(bry.funct.pca, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(bry.funct.pca, display="species", choices=c(1), scaling=2)*1.1,
     scores(bry.funct.pca, display="species", choices=c(2), scaling=2)*1.1,
     labels=rownames(scores(bry.funct.pca, display="species", scaling=2)),
     col="black", cex=1)   
dev.off()


############################################################################
####Figure S5. Functional trait frequencies along latitudinal gradients ####
############################################################################
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water) #select environmental data

my.formula <- y ~ x #formula used for plotting of trait frequencies as a function of latitude (linear)

################################
#### Vasc traits - latitude ####
################################
env<-merge(env, cwm_vascu, by = "row.names") #merge trait frequencies (CWM matrix) with environmental data for easy plotting. 

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Herbaceous','Shrub','Tree'), value.name = 'CWM_life')

g.1 <- ggplot(env.long, aes(x=Latitude, y=CWM_life, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Life form", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,0.8))
g.1


env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Height 1','Height 2','Height 3'), value.name = 'CWM_height')

g.2 <- ggplot(env.long, aes(x=Latitude, y=CWM_height, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Height", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,0.8))
g.2

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Seed 1','Seed 2','Seed 3'), value.name = 'CWM_Seed')

g.3<- ggplot(env.long, aes(x=Latitude, y=CWM_Seed, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Seed weight", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,0.8))
g.3

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('SLA 1','SLA 2','SLA 3'), value.name = 'CWM_SLA')

g.4 <- ggplot(env.long, aes(x=Latitude, y=CWM_SLA, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="SLA", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,0.8)) 
g.4

svglite("figure_S5_vasc_traits.svg")
grid.arrange(g.1, g.2, g.3, g.4,ncol=2)
dev.off()

################################
#### Moss traits - latitude ####
################################
#### Calculate CWM ####
cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_bryo) <- c('Acrocarpous','Pleurocarpous','Sphagnum','Tuft','Turf','Weft','Shoot 1','Shoot 2','Shoot 3','Dioicous','Monoicous','Seta 0','Seta 2','Seta 3','Seta 4','Diplolepidous alternate','Diplolepidous perfect','Haplolepideous','Nematodontus','No peristome','Spore 1','Spore 2','Spore 3','No Tomentum', 'Tomentum')

env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)
env<-merge(env, cwm_bryo, by = "row.names")

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Acrocarpous','Pleurocarpous','Sphagnum'), value.name = 'CWM_group')

g.1 <- ggplot(env.long, aes(x=Latitude, y=CWM_group, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Bryophyte group", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.1

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Turf','Tuft','Weft'), value.name = 'CWM_growth')

g.2 <- ggplot(env.long, aes(x=Latitude, y=CWM_growth, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Growth form", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.2

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Shoot 1','Shoot 2','Shoot 3'), value.name = 'CWM_shoot')

g.3<- ggplot(env.long, aes(x=Latitude, y=CWM_shoot, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Shoot length", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.3

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Seta 0','Seta 2','Seta 3', 'Seta 4'), value.name = 'CWM_seta')

g.4 <- ggplot(env.long, aes(x=Latitude, y=CWM_seta, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Seta length", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.4

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Haplolepideous','Nematodontus','No peristome','Diplolepidous alternate','Diplolepidous perfect'), value.name = 'CWM_peristome')

g.5 <- ggplot(env.long, aes(x=Latitude, y=CWM_peristome, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Peristome type", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.5

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Spore 1','Spore 2','Spore 3'), value.name = 'CWM_spore')

g.6 <- ggplot(env.long, aes(x=Latitude, y=CWM_spore, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Spore size", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.6

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('Monoicous','Dioicous'), value.name = 'CWM_sex')

g.7 <- ggplot(env.long, aes(x=Latitude, y=CWM_sex, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Sexual condition", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.7

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('No Tomentum','Tomentum'), value.name = 'Tomentum')

g.8 <- ggplot(env.long, aes(x=Latitude, y=Tomentum, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Tomentum", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.8

svglite("figure_S6_bry_traits.svg")
grid.arrange(g.1, g.2, g.3, g.4,g.5,g.6,g.7,g.8,ncol=2)
dev.off()
