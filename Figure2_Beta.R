#Load packages
require(effects)
require(nlme)
require(emmeans)
require(ggplot2)
require(dplyr)
require(car)
require(gridExtra)
require(svglite)
require(sjPlot)
require(adespatial)
require(vegan)
require(FD)
set_theme(base = theme_classic())  

#Open data
sites <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/sites_data.csv")
vascu <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/vascu_data.csv")
bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/bryo_data.csv")

vascu <- data.frame(vascu[,-1], row.names = vascu$X)
bryo <- data.frame(bryo[,-1], row.names = bryo$X)

##################################
#### Beta diversity taxonomic ####
##################################

# Computation using beta.div {adespatial} on 
# Hellinger-transformed species data
vas.beta <- beta.div(vascu, method = "hellinger", nperm = 9999)
summary(vas.beta)
vas.beta$beta  # SSTotal and BDTotal 

# Bryophyte
bry.beta <- beta.div(bryo, method = "hellinger", nperm = 9999)
summary(bry.beta)
bry.beta$beta  # SSTotal and BDTotal

sites$LCBD.bry.t <- bry.beta$LCBD
sites$LCBD.vas.t <- vas.beta$LCBD

#####LMEM 
env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, LCBD.vas.t, LCBD.bry.t)

env.scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)

############################
#### Vascular taxo LCBD ####
############################

result <- lme(log(LCBD.vas.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
summary(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

vif(result) 

#Plot
result <- lme(log(LCBD.vas.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

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

#################################
####Bryophyte taxonomic LCBD ####
#################################

result <- lme(log(LCBD.bry.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env.scale) 
summary(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

vif(result)

#Plot
result <- lme(log(LCBD.bry.t) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

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

###################################
#### Beta diversity functional ####
###################################

traits_vasc <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_vasc_data.csv")
traits_vasc <- as.data.frame(unclass(traits_vasc), stringsAsFactors = TRUE)
traits_vasc <- data.frame(traits_vasc, row.names = traits_vasc$Code)
traits_vasc <- traits_vasc[order(traits_vasc$Code),]
traits_vasc <- dplyr::select(traits_vasc, Port, Height, Seed, SLA)

vascu <-vascu[,order(colnames(vascu))]
vascu<- as.matrix(vascu)

cwm_vascu <- functcomp(traits_vasc, vascu, CWM.type = 'all')
colnames(cwm_vascu) <- c('Herbaceous', 'Shrub', 'Tree','Height 1', 'Height 2', 'Height 3', 'Seed 1', 'Seed 2', 'Seed 3', 'SLA 1', 'SLA 2', 'SLA 3')

traits_vasc_hel <- decostand(cwm_vascu, "hellinger")

traits_bryo <- read.csv("C:/Users/Utilisateur/Documents/MAITRISE/Scripts/Data_repository/traits_bry_data.csv")
traits_bryo <- as.data.frame(unclass(traits_bryo), stringsAsFactors = TRUE)
traits_bryo <- data.frame(traits_bryo, row.names = traits_bryo$ï..Code)
traits_bryo <- traits_bryo[order(traits_bryo$ï..Code),]
traits_bryo <- dplyr::select(traits_bryo, Life_strategy, Bryophyte.group, Life_form)
sapply(traits_bryo, class)

bryo <-bryo[,order(colnames(bryo))]
bryo<- as.matrix(bryo)

cwm_bryo <- functcomp(traits_bryo, bryo, CWM.type = 'all') # returns the frequencies of each class
colnames(cwm_bryo) <- c('Dominant', 'Perennial', 'Shuttle', 'Acrocarpous', 'Pleurocarpous', 'Sphagnum', 'Tuft', 'Turf', 'Weft')

traits_bryo_hel <- decostand(cwm_bryo, "hellinger")

# Computation using beta.div {adespatial} on 
# Hellinger-transformed species data
vas.beta.f <- beta.div(cwm_vascu, method = "hellinger", nperm = 9999)
summary(vas.beta.f)
vas.beta.f$beta  # SSTotal and BDTotal 

# Which species have a SCBD larger than the mean SCBD?
vas.beta.f$SCBD[vas.beta.f$SCBD >= mean(vas.beta.f$SCBD)]

# Bryophyte
bry.beta.f <- beta.div(cwm_bryo, method = "hellinger", nperm = 9999)
summary(bry.beta.f)
bry.beta.f$beta  # SSTotal and BDTotal

sites$LCBD.bry.funct <- bry.beta.f$LCBD
sites$LCBD.vas.funct <- vas.beta.f$LCBD

env <- dplyr::select(sites, SITES, Habitat, Latitude, Longitude, bioclim_1, bioclim_12, Thickness, Surface_water, LCBD.vas.funct, LCBD.bry.funct)

env_scale <- env
env.scale[,c(3:8)] <- scale(env.scale[,c(3:8)], center = TRUE, scale = TRUE)

#############################
#### Vasculat funct LCBD ####
#############################

result <- lme(log(LCBD.vas.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env_scale) 
summary(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

vif(result)

#Plot
result <- lme(log(LCBD.vas.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

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

##############################
#### Bryophyte funct LCBD ####
##############################

result <- lme(log(LCBD.bry.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env_scale) 
summary(result)

#Conditions
# normal plot of standardized residuals by habitat
qqnorm(result, ~ resid(., type = "p") | Habitat, abline = c(0, 1))
# normal plots of random effects
qqnorm(result, ~ranef(.))

vif(result)

#Plot
result <- lme(log(LCBD.bry.funct) ~ Habitat + Latitude + Longitude + bioclim_1 + bioclim_12 + Thickness + Surface_water + Latitude:Habitat, random = ~ 1|SITES, data = env)
CL <- emtrends(result, pairwise ~ Habitat, var="Latitude") 
test(CL)

df <- ggeffect(result, terms = c("Latitude","Habitat")) 
df$predicted <- exp(df$predicted)
df$conf.low <- exp(df$conf.low)
df$conf.high<-exp(df$conf.high)

mod_bry_lcbd_f <- ggplot(df, aes(x, predicted)) + 
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


svglite("figure2_beta_aleatoire.svg")
gridExtra::grid.arrange(mod_vas_lcbd_t, mod_bry_lcbd_t, mod_vas_lcbd_f,mod_bry_lcbd_f, ncol=2)
dev.off()

