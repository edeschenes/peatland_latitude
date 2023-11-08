
#Load packages
require(ggeffects)
require(nlme)
require(emmeans)
require(ggplot2)
require(dplyr)
require(car)
require(gridExtra)
require(svglite)
require(adespatial)
require(vegan)
require(FD)

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


#######################################################################
#### Calculate LCBD indices with site x species matrix (taxonomic) ####
#######################################################################

# Computation of LCBD using beta.div {adespatial} on Hellinger-transformed species data

#Vascular
vas.beta <- beta.div(vascu, method = "hellinger", nperm = 9999)
summary(vas.beta)
vas.beta$beta  

# Bryophyte
bry.beta <- beta.div(bryo, method = "hellinger", nperm = 9999)
summary(bry.beta)
bry.beta$beta  # SSTotal and BDTotal

#Extract LCBD indices and assign to each plot in the site x variables matrix
sites$LCBD.bry.t <- bry.beta$LCBD
sites$LCBD.vas.t <- vas.beta$LCBD

#Select only environmental variables used in the LCBD analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water, vascular richness, bryophyte richness)
env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum, LCBD.vas.t, LCBD.bry.t)

#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center). The unstandardized data will be used to obtain predicted values from the models to be used in the plots (Figure 2).
env.scale <- env
env.scale[,c(3:9)] <- scale(env.scale[,c(3:9)], center = TRUE, scale = TRUE)


##################################################################################
#### Linear mixed models - Vascular taxonomic LCBD as a function of variables ####
##################################################################################
#Log transform so conditions are met
result <- lme(log(LCBD.vas.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 

#Coefficients and Pvalues in Table 1 are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of LCBD as a function of latitude for bog/fens separately (used in Figure 2), use the non standardized data. This is also used for plotting of predicted values. 
result <- lme(log(LCBD.vas.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env)

#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately (values used in Figure 1)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Transform to exponential as the model was run on log transformed vascular richness
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

#Plot of vascular taxonomic LCBD as a function of latitude in bogs and fens (Figure 3A)
mod_vas_lcbd_t <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("solid", "solid"))+  
  labs(y="LCBD", x = "Latitude (degrees)")+
  ylim(0.001,0.0065)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

##############################################################################
#### Linear mixed models - Moss taxonomic LCBD as a function of variables ####
##############################################################################
#Log transform so conditions are met
result <- lme(log(LCBD.bry.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
#Coefficients and Pvalues in Table 1 are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of LCBD as a function of latitude for bog/fens separately (used in Figure 2), use the non standardized data. This is also used for plotting of predicted values. 
result <- lme(log(LCBD.bry.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env)

#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately (values used in Figure 1)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Transform to exponential as the model was run on log transformed vascular richness
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

#Plot of vascular taxonomic LCBD as a function of latitude in bogs and fens (Figure 3B)
mod_bry_lcbd_t <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("solid", "solid"))+  
  labs(y="LCBD", x = "Latitude (degrees)")+
  ylim(0.001,0.0065)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

#################################################################################
#### Calculate LCBD indices with site x trait class abundance matrix (functional)
#################################################################################

##########################
#### Vascular species ####
##########################

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

#######################
#### Calculate CWM ####
#######################
cwm_vascu <- functcomp(traits_vasc, vascu, CWM.type = 'all')# returns the frequencies of each class

#Modify trait names for plotting
colnames(cwm_vascu) <- c('Herbaceous', 'Shrub', 'Tree','Height 1', 'Height 2', 'Height 3', 'Seed 1', 'Seed 2', 'Seed 3', 'SLA 1', 'SLA 2', 'SLA 3')

######################
#### Moss species ####
######################
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

cwm_bryo <- subset(cwm_bryo, select = -c(No_Tomentum, Monoicous))

############################################################
#### Calculate LCBD indices with CWM matrix (functional)####
############################################################

# Computation using beta.div {adespatial} on Hellinger-transformed frequencies for each functionl trait class (CWM)

#Vascular
vas.beta.f <- beta.div(cwm_vascu, method = "hellinger", nperm = 9999)
summary(vas.beta.f)
vas.beta.f$beta  # SSTotal and BDTotal

# Bryophyte
bry.beta.f <- beta.div(cwm_bryo, method = "hellinger", nperm = 9999)
summary(bry.beta.f)
bry.beta.f$beta  # SSTotal and BDTotal

#Extract LCBD indices and assign to each plot in the site x variables matrix
sites$LCBD.bry.funct <- bry.beta.f$LCBD
sites$LCBD.vas.funct <- vas.beta.f$LCBD

#Select only environmental variables used in the LCBD analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water, vascular richness, bryophyte richness)
env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Substratum, LCBD.vas.funct, LCBD.bry.funct)


#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center). The unstandardized data will be used to obtain predicted values from the models to be used in the plots (Figure 2).
env.scale <- env
env.scale[,c(3:9)] <- scale(env.scale[,c(3:9)], center = TRUE, scale = TRUE)


##################################################################################
#### Linear mixed models - Vascular taxonomic LCBD as a function of variables ####
##################################################################################
#Log transform so conditions are met
result <- lme(log(LCBD.vas.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
#Coefficients and Pvalues in Table 1 are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of LCBD as a function of latitude for bog/fens separately (used in Figure 2), use the non standardized data. This is also used for plotting of predicted values.
result <- lme(log(LCBD.vas.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env)
#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately (values used in Figure 1)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Transform to exponential as the model was run on log transformed vascular richness
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

#Plot of vascular functional LCBD as a function of latitude in bogs and fens (Figure 3C)
mod_vas_lcbd_f <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("dashed", "dashed"))+  
  labs(y="LCBD", x = "Latitude (degrees)")+
  ylim(0,0.007)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

###############################################################################
#### Linear mixed models - Moss functional LCBD as a function of variables ####
###############################################################################
#Log transform so conditions are met
result <- lme(log(LCBD.bry.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
#Coefficients and Pvalues in Table 1 are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of LCBD as a function of latitude for bog/fens separately (used in Figure 2), use the non standardized data. This is also used for plotting of predicted values.
result <- lme(log(LCBD.bry.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Substratum + Latitude:Habitat, random = ~ 1|SITES, data = env)

#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately (values used in Figure 1)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Transform to exponential as the model was run on log transformed vascular richness
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

#Plot of moss functional LCBD as a function of latitude in bogs and fens (Figure 3D)
mod_bry_lcbd_f <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+  
  scale_linetype_manual(values = c("solid", "solid"))+  
  labs(y="LCBD", x = "Latitude (degrees)")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

#######################################################
####Download SVG file for modification in Inkscape ####
#######################################################
svglite("figure3_lcbd.svg", width = 6.6, fix_text_size = FALSE, pointsize =8)
gridExtra::grid.arrange(mod_vas_lcbd_t, mod_bry_lcbd_t, mod_vas_lcbd_f,mod_bry_lcbd_f, ncol=2)
dev.off()

