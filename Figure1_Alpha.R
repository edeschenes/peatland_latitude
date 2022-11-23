
#Load packages
require(ggeffects)
require(nlme)
require(emmeans) 
require(ggplot2)
require(dplyr)
require(car)
require(gridExtra)
require(svglite)
require(FD)
require(sjPlot)
require(effects)

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

#Select only environmental variables used in the alpha diversity analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water, vascular richness, bryophyte richness)
env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Vas_Ric_Tot, Bry_Ric_Tot)

#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center). The unstandardized data will be used to obtain predicted values from the models to be used in the plots (Figure 1). 
env.scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)

####################################################################
#### Vascular richness as a function of environmental variables ####
####################################################################

#Model of vascular richness as a function of standardized environmental variables, including a mixed effect for plots within sites
result <- lme(log(Vas_Ric_Tot) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale)

#Coefficients and Pvalues in Table 1 are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of richness as a function of latitude for bog/fens separately (used in Figure 1), use the non standardized data. This is also used for plotting of predicted values. 
result <- lme(log(Vas_Ric_Tot) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)

#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately (values used in Figure 1)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Transform to exponential as the model was run on log transformed vascular richness
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

#Plot of vascular richness as a function of latitude in bogs and fens (Figure 1A)
mod_vas_ric <- ggplot(df, aes(x, predicted)) + 
geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("dashed","solid"))+  
  labs(y="Richness", x = "Latitude (degrees)")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

################################################################
#### Moss richness as a function of environmental variables ####
################################################################

#Model of moss richness as a function of standardized environmental variables, including a mixed effect for sites
result <- lme(Bry_Ric_Tot ~ Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 

#Coefficients and Pvalues are obtained from the following:
summary(result)
anova(result)
#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of richness as a function of latitude for bog/fens separately, use the non standardized data in the model of moss richness as a function of environmental variables. This is used for plotting of predicted values. 
result <- lme(Bry_Ric_Tot ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)

#Determine the p-value and coefficient of richness as a function of latitude for bogs and fens separately.
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot of moss richness as a function of latitude in bogs and fens (Figure 1B)
mod_bry_ric <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("dashed", "dashed"))+  
  labs(y="Richness", x = "Latitude (degrees)")+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))


##########################################################################
####Calculation of functional alpha diversity : Functional dispersion ####
##########################################################################

require(FD)

#########################################
#### Calculate vascular species FDis ####
#########################################

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
vascu <- as.matrix(vascu %>% mutate_if(is.character,as.numeric))

#Calculate functional diversity indices for vascular functional diversity. Inputs are the species x trait matrix and the site x species matrix. 
fd_v <- dbFD(traits_vasc, vascu, corr = "cailliez", CWM.type = 'all',  ord = 'podani')

#Extract FDis from the calculation of functional diversity indices, and assign each FDis to it's plot in the site x environmental variables matrix
sites$FDis_v <- fd_v$FDis

#####################################
#### Calculate moss species FDis ####
#####################################

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
traits_bryo <- dplyr::select(traits_bryo, Bryophyte_group, Life_form, Shoot_length, Sexual_condition, Seta_length,Peristome_type, Spore_size, Tomentum)

#Order site x moss abundance matrix by column names (which are species codes). This is done to ensure that the species x trait matrix and the site x species matrix following the same order of species.
bryo <-bryo[,order(colnames(bryo))]

#Create site x vascular abundance matrix
bryo <- as.matrix(bryo %>% mutate_if(is.integer,as.numeric))

#Calculate functional diversity indices for moss functional diversity. Inputs are the species x trait matrix and the site x species matrix. 
fd_b <- dbFD(traits_bryo, bryo, corr = "cailliez", CWM.type = 'all', ord = 'podani')

#Extract FDis from the calculation of functional diversity indices, and assign each FDis to it's plot in the site x environmental variables matrix
sites$FDis_b <- fd_b$FDis

#Select only environmental variables used in the alpha diversity analyses (Habitat, latitude, longitude, mean annual temperature, mean annual precipitations, peat thickness, surface water, vascular richness, bryophyte richness) as well as FDis
env <- dplyr::select(sites, SITES,Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, FDis_v, FDis_b)

#Create a new data frame for standardized environmental variables. Standardize environmental variables (scale and center). The unstandardized data will be used to obtain predicted values from the models to be used in the plots (Figure 1). 
env.scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)


################################################################
#### Vascular FDis as a function of environmental variables ####
################################################################

#Model of vascular FDis as a function of standardized environmental variables, including a mixed effect for sites
result <- lme(FDis_v ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat,  random = ~ 1|SITES, data = env.scale) 

#Coefficients and Pvalues are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of richness as a function of latitude for bog/fens separately, use the non standardized data in the model of moss richness as a function of environmental variables. This is used for plotting of predicted values.
result <- lme(FDis_v ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat,  random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot of vascular FDis as a function of latitude in bogs and fens (Figure 1C)
mod_vas_fdis <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("dashed", "dashed"))+  
  labs(y="FDis", x = "Latitude (degrees)")+
  ylim(-0.1,0.5)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

############################################################
#### Moss FDis as a function of environmental variables ####
############################################################

result <- lme(FDis_b ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
#Coefficients and Pvalues are obtained from the following:
summary(result)

#VIFs were verified and always less than 20. 
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of richness as a function of latitude for bog/fens separately, use the non standardized data in the model of moss richness as a function of environmental variables. This is used for plotting of predicted values.
result <- lme(FDis_b ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env) 
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

#Obtain predicted values and confidence intervals for plotting purposes
df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot of moss richness as a function of latitude in bogs and fens (Figure 1D)
mod_bry_fdis <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("solid", "solid"))+  
  labs(y="FDis", x = "Latitude (degrees)")+
  ylim(-0.1,0.5)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

#######################################################
####Download SVG file for modification in Inkscape ####
#######################################################
svglite("figure1_alpha_aleatoire.svg")
grid.arrange(mod_vas_ric, mod_bry_ric, mod_vas_fdis, mod_bry_fdis,ncol=2)
dev.off()

