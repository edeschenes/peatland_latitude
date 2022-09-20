
#Load packages
require(eulerr)
require(vegan) 
require(FD)
require(svglite)

#Open data
sites <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/sites_data.csv")
vascu <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/vascu_data.csv")
bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/bryo_data.csv")

vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)

###################
####Composition####
###################
vas.hel <- (decostand(vascu, "hellinger"))
bry.hel <- (decostand(bryo, "hellinger"))

sites <- data.frame(sites, row.names = sites$ID_Quadrat)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

env[,c(2:7)] <- scale(env[,c(2:7)], center = TRUE, scale = TRUE)
env$Habitat <- factor(env$Habitat)

####Vascular####
vas.rda <- rda(vas.hel ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, env) 
(vas.rda <- rda(vas.hel ~ ., env))
summary(vas.rda)	# Scaling 2 (default) 

plot(vas.rda) # use base plot, might be done with ggplot2
anova(vas.rda) # is the model significant?


####Bryophyte####
bry.rda <- rda(bry.hel ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, env) 
(bry.rda <- rda(bry.hel ~ ., env)) 
summary(bry.rda)	# Scaling 2 (default) 

plot(bry.rda) # use base plot, might be done with ggplot2
anova(bry.rda) # is the model significant?


####Variation partitioning####

bioclim <- env[, c(4:5)]
space <- env[, c(2:3)]
local <- env[, c(1,6,7)]

(spe.part.vas <- varpart(vascu, bioclim,space, local, transfo= 'hel'))
plot(spe.part.vas, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim', 'Space', 'Local'))

venn_vasc_tax <- euler(c(A = 1,          # Draw pairwise venn diagram
                     B = 1.4,
                     C = 12.4,
                     "A&B" = 2.1,
                     'A&C'= 0.5,
                     'B&C' = 0.5,
                     'A&B&C' = 0))
p_venn_vasc_tax <- plot(venn_vasc_tax, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)
p_venn_vasc_tax

(frac.1 <- anova(rda(vas.hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(vas.hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(vas.hel, local, cbind.data.frame(bioclim, space)), step=10000))


(spe.part.bry <- varpart(bryo, bioclim,space, local,transfo= 'hel'))
plot(spe.part.bry, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim', 'Space', 'Local'))

venn_bryo_tax <- euler(c(A = 0.6,          # Draw pairwise venn diagram
                         B = 0.9,
                         C = 6.9,
                         "A&B" = 1.9,
                         'A&C'= 0.3,
                         'B&C' = 0.3,
                         'A&B&C' = 0.5))
p_venn_bryo_tax <- plot(venn_bryo_tax, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)
p_venn_bryo_tax

(frac.1 <- anova(rda(bry.hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(bry.hel, local, cbind.data.frame(bioclim, space)), step=10000))
(frac.3 <- anova(rda(bry.hel, space, cbind.data.frame(bioclim, local)), step=10000))

##################
####Functional####
##################

traits_vasc <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_vasc_data.csv")
traits_vasc <- as.data.frame(unclass(traits_vasc), stringsAsFactors = TRUE)
traits_vasc <- data.frame(traits_vasc, row.names = traits_vasc$Code)
traits_vasc <- traits_vasc[order(traits_vasc$Code),]
traits_vasc <- dplyr::select(traits_vasc, Port, Height, Seed, SLA)

vascu <-vascu[,order(colnames(vascu))]
vascu<- as.matrix(vascu)

cwm_vascu <- functcomp(traits_vasc, vascu, CWM.type = 'all')
colnames(cwm_vascu)
colnames(cwm_vascu) <- c('Herbaceous', 'Shrub', 'Tree','Height 1', 'Height 2', 'Height 3', 'Seed 1', 'Seed 2', 'Seed 3', 'SLA 1', 'SLA 2', 'SLA 3')

traits_vasc_hel <- decostand(cwm_vascu, "hellinger")

traits_bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_bry_data.csv")
traits_bryo <- as.data.frame(unclass(traits_bryo), stringsAsFactors = TRUE)
traits_bryo <- data.frame(traits_bryo, row.names = traits_bryo$Code)
traits_bryo <- traits_bryo[order(traits_bryo$Code),]
traits_bryo <- dplyr::select(traits_bryo, Life_strategy, Bryophyte.group, Life_form)
sapply(traits_bryo, class)


bryo <-bryo[,order(colnames(bryo))]
bryo<- as.matrix(bryo)

cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class
colnames(cwm_bryo)
colnames(cwm_bryo) <- c('Dominant', 'Perennial', 'Shuttle',  'Acrocarpous', 'Pleurocarpous', 'Sphagnum', 'Tuft', 'Turf', 'Weft')

traits_bryo_hel <- decostand(cwm_bryo, "hellinger")


####Vascular####
(vas.rda <- rda(traits_vasc_hel ~ ., env))
summary(vas.rda) # Scaling 2 (default) 

plot(vas.rda) # use base plot, might be done with ggplot2
anova(vas.rda) # is the model significant?

(spe.part.vas <- varpart(traits_vasc_hel, bioclim,space, local))
plot(spe.part.vas, digits = 2, bg = c("red", "blue", 'green','yellow'),Xnames = c('Bioclimatique', 'Space', 'Local'))

(frac.1 <- anova(rda(traits_vasc_hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(traits_vasc_hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(traits_vasc_hel, local, cbind.data.frame(bioclim, space)), step=10000))


venn_vasc_funct <- euler(c(A = 0.5,          # Draw pairwise venn diagram
                           B = 1.6,
                           C = 13.8,
                           "A&B" = 2.7,
                           'A&C'= 0.1,
                           'B&C' = 0.5,
                           'A&B&C' =0))
p_venn_vasc_funct<- plot(venn_vasc_funct, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)
p_venn_vasc_funct


####Bryophyte####
bry.rda <- rda(traits_bryo_hel ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, env) 
(bry.rda <- rda(traits_bryo_hel ~ ., env)) 
summary(bry.rda)	# Scaling 2 (default) 

plot(bry.rda) # use base plot, might be done with ggplot2
anova(bry.rda) # is the model significant?

(spe.part.bry <- varpart(traits_bryo_hel, bioclim,space, local))
plot(spe.part.bry, digits = 2, bg = c("red", "blue", 'green'),Xnames = c('Bioclim','Space', 'Local'))

(frac.1 <- anova(rda(traits_bryo_hel, bioclim, cbind.data.frame(space, local)), step=10000))
(frac.2 <- anova(rda(traits_bryo_hel, space, cbind.data.frame(bioclim, local)), step=10000))
(frac.3 <- anova(rda(traits_bryo_hel, local, cbind.data.frame(bioclim, space)), step=10000))


venn_bryo_funct <- euler(c(A = 1,          # Draw pairwise venn diagram
                           B = 1.1,
                           C = 12.8,
                           "A&B" = 4.7,
                           'A&C'= 1.1,
                           'B&C' = 0.8,
                           'A&B&C' = 0))
p_venn_bryo_funct <- plot(venn_bryo_funct, quantities = list(type = c("counts"), cex = seq(1, 1, length.out = 4)), labels = list(labels = c("Bioclimatic", 'Space', 'Local'),  cex = seq(1,1, length.out = 4)),adjust_labels =TRUE)
p_venn_bryo_funct


gridExtra::grid.arrange(p_venn_vasc_funct, p_venn_bryo_funct, ncol=2)

svglite("figure3_varpart.svg")
gridExtra::grid.arrange(p_venn_vasc_tax, p_venn_bryo_tax,p_venn_vasc_funct, p_venn_bryo_funct, ncol=2)
dev.off()


####################
#### RDA TRAITS ####
####################

sites <- data.frame(sites, row.names = sites$ID_Quadrat)
env <- dplyr::select(sites, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water)

env[,c(2:7)] <- scale(env[,c(2:7)], center = TRUE, scale = TRUE)
colnames(env) <- c("Habitat", "Latitude",'Longitude', 'Annual temperature', 'Annual precipitations', 'Thickness', 'Surface water')

env$Habitat=as.factor(env$Habitat)
par(mfrow=c(1,1))
scl <- 2 ## scaling == 3
colvec <- c('lightcoral','lightblue3')
pointvec <- c(19,17)


####Vascular####
vas.rda <-rda(traits_vasc_hel~ ., env) # Observe the shortcut formula
summary(vas.rda)

svglite("figure4_RDA_vasc.svg")
plot(vas.rda, scaling=2, main="Triplot RDA - vascular traits", type="none", xlab=c("RDA1 (89.53%)"), ylab=c("RDA2 (4.58%)"), xlim=c(-1, 1), ylim=c(-1.4,1.4), cex =1)
with(env, points(vas.rda, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec, pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(vas.rda, display="species", choices=c(1), scaling=2)*2,
       scores(vas.rda, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(vas.rda, display="species", choices=c(1), scaling=2)*2.4,
     scores(vas.rda, display="species", choices=c(2), scaling=2)*2.4,
     labels=rownames(scores(vas.rda, display="species", scaling=2)),
     col="black", cex=1)   
arrows(0,0,
       scores(vas.rda, display="bp", choices=c(1), scaling=2)*2,
       scores(vas.rda, display="bp", choices=c(2), scaling=2)*2,
       col="red",length = 0.1)
text(scores(vas.rda, display="bp", choices=c(1), scaling=2)*2,
     scores(vas.rda, display="bp", choices=c(2), scaling=2)*2,
     labels=rownames(scores(vas.rda, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1) 
with(env,
     points(scores(vas.rda, display = "cn", pch = 10)*2))
with(env,
    points(scores(vas.rda, display = "cn", pch = 10)*2))

dev.off()


####Bryophyte####

bry.rda <-rda(traits_bryo_hel~ ., env) # Observe the shortcut formula
summary(bry.rda)	# Scaling 2 (default) 

svglite("figure4_RDA_bryo.svg")
plot(bry.rda, scaling=2, main="Triplot RDA - bryo traits", type="none", xlab=c("RDA1 (57.61%)"), ylab=c("RDA2 (35.33%)"),xlim=c(-1.5,1.5), ylim=c(-0.75,1.1), cex = 1)
with(env, points(bry.rda, display = "sites", col = colvec[Habitat],
                 scaling = scl, pch = pointvec[Habitat], bg = colvec[Habitat]))
with(env, legend("topright", legend = levels(Habitat), bty = "n",
                 col = colvec,pch = pointvec[Habitat], pt.bg = colvec))
arrows(0,0,
       scores(bry.rda, display="species", choices=c(1), scaling=2)*2,
       scores(bry.rda, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(bry.rda, display="species", choices=c(1), scaling=2)*2.1,
     scores(bry.rda, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(bry.rda, display="species", scaling=2)),
     col="black", cex=1)   
arrows(0,0,
       scores(bry.rda, display="bp", choices=c(1), scaling=2)*1.5,
       scores(bry.rda, display="bp", choices=c(2), scaling=2)*1.5,
       col="red",length = 0.1)
text(scores(bry.rda, display="bp", choices=c(1), scaling=2)*1.5+0.05,
     scores(bry.rda, display="bp", choices=c(2), scaling=2)*1.5+0.05,
     labels=rownames(scores(bry.rda, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)
text(scores(bry.rda, display="bp", choices=c(1), scaling=2)+0.05,
     scores(bry.rda, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(bry.rda, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1) 
dev.off()


