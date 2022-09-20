
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
set_theme(base = theme_classic())  

#Open data
sites <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/sites_data.csv")
vascu <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/vascu_data.csv")
bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/bryo_data.csv")

vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)

sites <- data.frame(sites, row.names = sites$ID_Quadrat)
env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, Vas_Ric_Tot, Bry_Ric_Tot)

env.scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)

###########################
#### VASCULAR RICHNESS ####
###########################

result <- lme(log(Vas_Ric_Tot) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale)
summary(result)
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of bog/fens seperately, use the non standardized data 
result <- lme(log(Vas_Ric_Tot) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

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

############################
#### BRYOPHYTE RICHNESS ####
############################

result <- lme(Bry_Ric_Tot ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
summary(result)
anova(result)
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of bog/fens seperately, use the non standardized data 
result <- lme(Bry_Ric_Tot ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot
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

#########################
####Functional alpha ####
#########################

require(FD)

####Vascular #### 
traits_vasc <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_vasc_data.csv")
traits_vasc <- as.data.frame(unclass(traits_vasc), stringsAsFactors = TRUE)
traits_vasc <- data.frame(traits_vasc, row.names = traits_vasc$Code)
traits_vasc <- traits_vasc[order(traits_vasc$Code),]
traits_vasc <- dplyr::select(traits_vasc, Port, Height, Seed, SLA)

vascu <-vascu[,order(colnames(vascu))]
vascu<- as.matrix(vascu)

fd_v <- dbFD(traits_vasc, vascu, corr = "cailliez", CWM.type = 'all')
sites$FDis_v <- fd_v$FDis

#### Bryophyte ####
traits_bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_bry_data.csv")
traits_bryo <- as.data.frame(unclass(traits_bryo), stringsAsFactors = TRUE)
traits_bryo <- data.frame(traits_bryo, row.names = traits_bryo$Code)
traits_bryo <- traits_bryo[order(traits_bryo$Code),]
traits_bryo <- dplyr::select(traits_bryo, Life_strategy, Bryophyte.group, Life_form)
sapply(traits_bryo, class)

bryo <-bryo[,order(colnames(bryo))]
bryo<- as.matrix(bryo)

fd_b <- dbFD(traits_bryo, bryo, corr = "cailliez", CWM.type = 'all')
sites$FDis_b <- fd_b$FDis

env <- dplyr::select(sites, SITES,Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, FDis_v, FDis_b)

env.scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)

#######################
#### Vascular FDis ####
#######################

result <- lme(FDis_v ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat,  random = ~ 1|SITES, data = env.scale) 
summary(result)
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of bog/fens seperately, use the non standardized data 
result <- lme(FDis_v ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat,  random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot
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

########################
#### Bryophyte FDis ####
########################

result <- lme(FDis_b ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
summary(result)
vif(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

#To calculate the p-value and coefficients of bog/fens seperately, use the non standardized data 
result <- lme(FDis_b ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env) 
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 

#Plot
mod_bry_fdis <- ggplot(df, aes(x, predicted)) + 
  geom_line(aes(linetype=group, color=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = c( "lightcoral","lightblue3"))+
  scale_fill_manual(values = c("lightcoral","lightblue3"))+ 
  scale_linetype_manual(values = c("dashed", "solid"))+  
  labs(y="FDis", x = "Latitude (degrees)")+
  ylim(-0.1,0.5)+
  theme(axis.title = element_text(size = 12),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12))

######################
####Final figures ####
######################
svglite("figure1_alpha_aleatoire.svg")
grid.arrange(mod_vas_ric, mod_bry_ric, mod_vas_fdis, mod_bry_fdis,ncol=2)
dev.off()

