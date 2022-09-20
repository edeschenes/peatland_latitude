#Scripts to create figures S2, S3 and S4. 
library(ggcorrplot)
require(GGally)

#Open data
sites <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/sites_data.csv")
vascu <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/vascu_data.csv")
bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/bryo_data.csv")

vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)

sites <- data.frame(sites, row.names = sites$ID_Quadrat)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)
env$Habitat=as.factor(env$Habitat)

###BOG
bog <- sites[sites[, 'Habitat'] == 'Bog', ] 
quad_bog <- unique(bog$ID_Quadrat)

bog_bry <- bryo[row.names(bryo) %in% quad_bog,]
bog_vas <- vascu[row.names(vascu) %in% quad_bog,]

#Remove species that are not present in bogs
bog_vas <- bog_vas[,colSums(bog_vas != 0) > 0]
bog_bry <- bog_bry[,colSums(bog_bry != 0) > 0]

env_bog<- env[row.names(env) %in% quad_bog,] 
env_bog[,c(2:7)] <- as.data.frame(scale(env_bog[,c(2:7)], center = TRUE, scale = TRUE))

###FEN
fen <- sites[sites[, 'Habitat'] == 'Fen', ] 
quad_fen <- unique(fen$ID_Quadrat)

fen_bry <- bryo[row.names(bryo) %in% quad_fen,]

fen_vas <- vascu[row.names(vascu) %in% quad_fen,]

env_fen<- env[row.names(env) %in% quad_fen,]
env_fen[,c(2:7)] <- as.data.frame(scale(env_fen[,c(2:7)], center = TRUE, scale = TRUE))

#Remove species that are not present in bogs
bog_vas <- bog_vas[,colSums(bog_vas != 0) > 0]
bog_bry <- bog_bry[,colSums(bog_bry != 0) > 0]

##################################################
#### Figure S2 - PCA of bioclimatic variables ####
##################################################

env1 <- dplyr::select(sites,bioclim_1, bioclim_2, bioclim_3, bioclim_4, bioclim_5, bioclim_6, bioclim_7, bioclim_8, bioclim_9, bioclim_10, bioclim_11, bioclim_12, bioclim_13, bioclim_14, bioclim_15, bioclim_16, bioclim_17, bioclim_18, bioclim_19)

env.pca <- rda(env1, scale = TRUE)
env.pca
summary(env.pca) # Default scaling 2

#Plot
par(mfrow=c(1,1))
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)
sites$Habitat=as.factor(sites$Habitat)
  
svglite("pca_bioclim.svg")
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
#### Figure S3 - Correlation plot of environmental variables ####
#################################################################

env_bog<- dplyr::select(env_bog, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)
colnames(env_bog) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water')

env_fen<- dplyr::select(env_fen, Latitude,Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)
colnames(env_fen) <- c("Latitude",'Longitude', 'Annual mean temperature', 'Annual precipitation','Peat thickness', 'Surface water')

r_bog <- cor(env_bog, use="complete.obs")
r_fen <- cor(env_fen, use="complete.obs")

cor_bog <- ggcorr(env_bog, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)
cor_fen <- ggcorr(env_fen, method = c("pairwise", "pearson"), label = TRUE, hjust = 0.75)

svglite("figureS3_corplot.svg")
gridExtra::grid.arrange(cor_bog, cor_fen, ncol=2)
dev.off()

##############################################################
#### Figure S4 - Environmental variables in bogs and fens ####
##############################################################
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

p1<-ggplot(data= env, aes(Habitat, Latitude))+
  stat_boxplot( aes(Habitat, Latitude), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Latitude),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Latitude (degrees)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Latitude~Habitat, data=env)
anova(fm1)

p2 <-ggplot(data= env, aes(Habitat, Longitude))+
  stat_boxplot( aes(Habitat, Longitude), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Longitude),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Longitude (degrees)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Longitude~Habitat, data=env)
anova(fm1)

p3 <-ggplot(data= env, aes(Habitat, bioclim_1))+
  stat_boxplot( aes(Habitat, bioclim_1), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, bioclim_1),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Annual temperature (°C)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(bioclim_1~Habitat, data=env)
anova(fm1)

p4<-ggplot(data= env, aes(Habitat, bioclim_12))+
  stat_boxplot( aes(Habitat, bioclim_12), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, bioclim_12),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Annual precipitation (mm)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(bioclim_12~Habitat, data=env)
anova(fm1)

p5<-ggplot(data= env, aes(Habitat, Thickness))+
  stat_boxplot( aes(Habitat, Thickness), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Thickness),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3) +
  labs(y="Peat thickness (cm)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Thickness~Habitat, data=env)
anova(fm1)

p6<-ggplot(data= env, aes(Habitat, Surface_water))+
  stat_boxplot( aes(Habitat, Surface_water), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Habitat, Surface_water),outlier.shape=1) +    
  stat_summary(fun=mean, geom="point", size=3)+
  labs(y="Surface water (%)", x = "Habitat")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

fm1 <- aov(Surface_water~Habitat, data=env)
anova(fm1)

svglite("figureS4_env.svg", width = 14, height = 7)
grid.arrange(p1,p3,p5,p2,p4,p6,ncol=3)
dev.off()

