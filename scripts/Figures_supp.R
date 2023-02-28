#Scripts to create supplementary figures 

#Load packages
library(ggcorrplot)
require(GGally)
require(vegan) 
require(FD)
require(svglite)
require(gridExtra)
library(reshape2)
require(ggpmisc)

set_theme(base = theme_classic()) 

#Load data
#Site x environmental variables matrix
sites <- read.csv("Data_repository/sites_data.csv")
#Site x species abundance matrix for vascular species
vascu <- read.csv("Data_repository/vascu_data.csv")
#Site x species abundance matrix for moss species
bryo <- read.csv("Data_repository/bryo_data.csv")

#Assign plotID to row names and remove plotID as a column
vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)
sites <- data.frame(sites, row.names = sites$ID_Quadrat)

#Select only environmental variables used in the supplementary analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum)

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
env_bog[,c(2:8)] <- as.data.frame(scale(env_bog[,c(2:8)], center = TRUE, scale = TRUE))

####Fens ####
fen <- sites[sites[, 'Habitat'] == 'Fen', ] #select fens
quad_fen <- unique(fen$ID_Quadrat) #assign fen quadrat IDs to object

fen_vas <- vascu[row.names(vascu) %in% quad_fen,] #select fens in site x vascular abundance matrix
fen_bry <- bryo[row.names(bryo) %in% quad_fen,] #select fens in site x moss abundance matrix

#Create environmental dataframe with only bogs
env_fen<- env[row.names(env) %in% quad_fen,]#select fens in site x enviro matrix
#Scale and center environmental variables
env_fen[,c(2:8)] <- as.data.frame(scale(env_fen[,c(2:8)], center = TRUE, scale = TRUE))


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
svglite("figure_S1_pca_bioclim.svg", width = 6.6, fix_text_size = FALSE, pointsize =8)
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


####################################
#### Figure S2 - PCOA BOG / FEN ####
####################################
#Betadisper: multivariate dispersion analysis
#permutest.betadisper: comparison of group mean dispersions
#rda and anova.cca: tests of centroid locations
bry.bray <- vegdist(bryo)
vas.bray <- vegdist(vascu)

gr <- env$Habitat

b.disp.bry <- betadisper(bry.bray, group = gr, type = c("centroid"), bias.adjust = FALSE)
b.disp.vas <- betadisper(vas.bray, group = gr, type = c("centroid"), bias.adjust = FALSE)

## S3 method for class 'betadisper'
anova(b.disp.bry) 
anova(b.disp.vas) 

## Permutation test for F
permutest(b.disp.bry, pairwise = TRUE, permutations = 9999) #Permdisp: p = 0.01
permutest(b.disp.vas, pairwise = TRUE, permutations = 9999) #Permdisp: p = 0.01

####Test of centroid locations 
adonis2(bry.bray ~ Habitat, data=env, permutations=9999) 
adonis2(vas.bray ~ Habitat, data=env, permutations=9999)

#Download PCoA biplots 
svglite("figureS2_pcoa.svg",  width = 6.6, fix_text_size = FALSE, pointsize =8)
par(mfrow=c(1,2))
plot(b.disp.vas, ellipse = TRUE, hull = FALSE, col = c('lightcoral','lightblue3')) 
plot(b.disp.bry, ellipse = TRUE, hull = FALSE,col = c('lightcoral','lightblue3')) 
dev.off()

#Distance to centroid figures. Added manually in inkscape to the PCoA plots
svglite("figureS2_distcent.svg")
par(mfrow=c(1,2))
boxplot(b.disp.vas, ylab = "Distance to centroid", main = 'Vascular',col = c("lightcoral", "lightblue3","lightcoral", "lightblue3","lightcoral", "lightblue3"), cex.axis =1.5 ,cex.lab=1.5)
boxplot(b.disp.bry, ylab = "Distance to centroid", main = 'Bryophyte',col =  c("lightcoral", "lightblue3","lightcoral", "lightblue3","lightcoral", "lightblue3"), cex.axis = 1.5,cex.lab=1.5)
dev.off()

################################################
#### Figure S3 - PCA of species composition ####
################################################
#Run PCA of species
vas.hel <- (decostand(vascu, "hellinger"))
bry.hel <- (decostand(bryo, "hellinger"))

vas.pca <- rda(vas.hel) 
summary(vas.pca)	# Scaling 2 (default) 

bry.pca <- rda(bry.hel) 
summary(bry.pca)	# Scaling 2 (default) 

#Prep for PCA plot
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)
sites$Habitat=as.factor(sites$Habitat)

svglite("figure_S3_pca_species.svg",  width = 6.6, fix_text_size = FALSE, pointsize =8)
par(mfrow=c(2,1))
plot(vas.pca, scaling=2,, type="none", xlab=c("PCA1 (22.25%)"), ylab=c("PCA2 (11.51%)"))
with(env, points(vas.pca, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(vas.pca, display="species", choices=c(1), scaling=2),
       scores(vas.pca, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(vas.pca, display="species", choices=c(1), scaling=2),
     scores(vas.pca, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(vas.pca, display="species", scaling=2)),
     col="black", cex=1)  

plot(bry.pca, scaling=2, type="none", xlab=c("PCA1 (18.16%)"), ylab=c("PCA2 (10.80%)"))
with(env, points(bry.pca, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(bry.pca, display="species", choices=c(1), scaling=2),
       scores(bry.pca, display="species", choices=c(2), scaling=2),
       col="black",length=0)
text(scores(bry.pca, display="species", choices=c(1), scaling=2),
     scores(bry.pca, display="species", choices=c(2), scaling=2),
     labels=rownames(scores(bry.pca, display="species", scaling=2)),
     col="black", cex=1)   
dev.off()

#################################################################
#### Figure S4 - Correlation plot of environmental variables ####
#################################################################

#Select environmental variables in bogs
env_bog<- dplyr::select(env_bog, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum)

#Change column names for plotting
colnames(env_bog) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water', 'Substratum')

#Select environmental variables in fens
env_fen<- dplyr::select(env_fen, Latitude,Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum)

#Change column names for plotting
colnames(env_fen) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water', 'Substratum')

#Calculate Pearson correlations between each pair of environmental variables
r_bog <- cor(env_bog, use="complete.obs")
r_fen <- cor(env_fen, use="complete.obs")

#Plot pairwise plot of Pearson correlations between environmental variables 
cor_bog <- ggcorr(env_bog, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)
cor_fen <- ggcorr(env_fen, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)

####Download SVG file for modification in Inkscape ####
svglite("figure_S2_corplot.svg", width = 6.6, fix_text_size = FALSE, pointsize =8)
gridExtra::grid.arrange(cor_bog, cor_fen, ncol=2)
dev.off()


##############################################################
#### Figure S5 - Environmental variables in bogs and fens ####
##############################################################
#Select environmental variables
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum)

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
anova(fm1) #P-values were used to determine significant differences between bogs and fens, and letters were added in Inkscape.

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
anova(fm1) #P-values were used to determine significant differences between bogs and fens, and letters were added in Inkscape.

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
anova(fm1) #P-values were used to determine significant differences between bogs and fens, and letters were added in Inkscape.

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
anova(fm1) #P-values were used to determine significant differences between bogs and fens, and letters were added in Inkscape.

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

p7<-ggplot(data= env, aes(Habitat, Substratum))+
  stat_boxplot( aes(Habitat, Deposit), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Substratum),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3)+
  labs(y="Substratum", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Substratum~Habitat, data=env)#run anova to determine difference between bogs and fens. 
anova(fm1) #P-values were used to determine significatn differences between bogs and fens, and letters were added in Inkscape.

####Download SVG file for modification in Inkscape ####
svglite("figure_S3_env.svg", width = 6, height = 12)
grid.arrange(p1,p3,p5,p2,p4,p6, p7, ncol=2)
dev.off()

############################################
####Figure S6. PCA of functional traits ####
############################################
#PCA is conducted on the site x trait class frequencies matrix.
#Step 1. Calculate trait class frequencies with function functcomp
#Step 2. PCA on Hellinger transformed site x trait class frequencies matrix

#Open vascular trait data 
traits_vasc <- read.csv("Data_repository/traits_vasc_data.csv",sep = ";", header= TRUE)

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
colnames(cwm_vascu) <- c('Herbaceous', 'Shrub', 'Tree','Height 1', 'Height 2', 'Height 3','Seed 1', 'Seed 2', 'Seed 3', 'SLA 1', 'SLA 2', 'SLA 3')

####################################
#### Calculate moss species CWM ####
####################################
#Open moss trait data 
traits_bryo <- read.csv("Data_repository/traits_bry_data.csv", sep = ",", header= TRUE)

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
traits_bryo <- dplyr::select(traits_bryo, Bryophyte.group, Life_form, Shoot_length, Sexual_condition, Seta_length,Peristome_type, Spore_size, Tomentum)

#Order site x moss abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species.
bryo <-bryo[,order(colnames(bryo))]

#Create site x vascular abundance matrix
bryo<- as.matrix(bryo)

#######################
#### Calculate CWM ####
#######################

cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_bryo) <- c('Acrocarpous','Pleurocarpous','Sphagnum','Tuft','Turf','Weft','Shoot 1','Shoot 2','Shoot 3','Dioicous','Monoicous','Seta 0','Seta 2','Seta 3','Seta 4','No peristome','Peristome perfect','Peristome specialized','Spore 1','Spore 2','Spore 3','No_Tomentum', 'Tomentum')

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

traits_bryo_hel <- sqrt(cwm_bryo)#Hellinger transformation on site x trait class frequencies
bry.funct.pca <- rda(traits_bryo_hel) #perform PCA for moss traits
summary(bry.funct.pca)

#Plot
svglite("figure_S6_pca_traits.svg", width= 6.6, fix_text_size = FALSE, pointsize =8.5)
par(mfrow = c(2,1))
plot(vas.funct.pca, scaling=2, type="none", xlab=c("PCA1 (56.02%)"), ylab=c("PCA2 (12.92%)"), xlim=c(-1.3, 1.3), ylim=c(-1,1), cex =1)
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
plot(bry.funct.pca, scaling=2,type="none", xlab=c("PCA1 (41.26%)"), ylab=c("PCA2 (19.92%)"), xlim=c(-1.3, 1.3), ylim=c(-1.5,1.5), cex = 1)
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
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum) #select environmental data

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

####Download SVG file for modification in Inkscape ####
svglite("figure_S7_vasc_traits.svg", , width= 6.6, fix_text_size = FALSE, pointsize =8.5)
grid.arrange(g.1, g.2, g.3, g.4,ncol=2)
dev.off()

################################
#### Moss traits - latitude ####
################################

env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum)
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

###stat_poly_eq not working for bryophyte.group. We manually calculated the values using lm, since that is what stat_poly_eq does. 
acrocarpous <- env.long[env.long$variable == 'Acrocarpous', ]
pleurocarpous <- env.long[env.long$variable == 'Pleurocarpous', ]
sphagnum <- env.long[env.long$variable == 'Sphagnum', ]

summary(lm(CWM_group ~ Latitude, data = acrocarpous))
summary(lm(CWM_group ~ Latitude, data = pleurocarpous))
summary(lm(CWM_group ~ Latitude, data = sphagnum))
#These values were manually added to the plot in inkscape.
####

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

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('No peristome','Peristome perfect','Peristome specialized'), value.name = 'CWM_peristome')

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

env.long <- melt(env, id = c("Row.names", 'Latitude'), measure = c('No_Tomentum','Tomentum'), value.name = 'Tomentum')

g.8 <- ggplot(env.long, aes(x=Latitude, y=Tomentum, color=variable)) +
  geom_smooth(method=lm, se=TRUE, formula = my.formula, aes(linetype=variable, color=variable)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~~~")), parse = TRUE, size = 3)+
  labs(title = '',y="Tomentum", x = "Latitude (degrees)")+
  theme_classic(base_size = 10)+
  theme(legend.position =  c(0.9, 0.9))+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(0,1.1))
g.8


####Download SVG file for modification in Inkscape ####
svglite("figure_S8_bry_traits.svg", width= 6.6,height = 9, fix_text_size = FALSE, pointsize =8.5)
grid.arrange(g.1, g.2, g.3, g.4,g.5,g.6,g.7,g.8,ncol=2)
dev.off()
