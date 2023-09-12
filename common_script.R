# Project 7 - How is invertebrate composition affected by elevation in Glenmore?
# Field course - September 2023
# Common script
# Adela, Adela, Emily, Lubomir, Nina

library(tidyverse)

meall <- read.csv("meall_data.csv")
meall <- meall[1:20,]

## Question 1 - inv affected by elevation -------------------------------

# richness
# As showed by models in question 3, richness is not affected by elevation in linear models.
par(mfrow = c(1,1))
plot(inv_order_rich ~ elev_m, data = meall)
# This graph shows that possibly, higher richness can be found in mid altitudes.
# We can't test that with a linear model, so might go beyond the scope of our stat skills.

# diversity
# As in question 3, diversity was not explained by elevation.
plot(shannon_diversity ~ elev_m, data = meall)  # nothing
# However, it seems to be well explained by vegetation, particularly by moss cover.
plot(shannon_diversity ~ moss_cov, data = meall)

# abundance
# As in question 3, elevation explains some variance in abundance, however, moss cover
# is a better predictor.
plot(inv_abun ~ elev_m, data = meall)
plot(inv_abun ~ moss_cov, data = meall)

## Question 2 - plants affected by elevation ------------------------------

# plant richness ~ elevation
ggplot(bada, aes(x = elev_m, y = plant_sp_rich)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()


# heather ~ elevation
ggplot(bada[2:20,], aes(x = elev_m, y = heather_cov)) +  # plot 1 is outlier, so exluded
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

# moss ~ elevation
ggplot(bada, aes(x = elev_m, y = moss_cov)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

par(mfrow = c(2,2))
lm <- lm(log(moss_cov+1) ~ elev_m, data = bada)
summary(lm)
lm0 <- lm(log(moss_cov+1) ~ 1, data = bada)
anova(lm, lm0, test = 'Chi')
plot(lm)
par(mfrow=c(1,1))
library(DHARMa)
sims <- simulateResiduals(lm)
plot(sims)

# lichen ~ elevation
ggplot(bada, aes(x = elev_m, y = lichen_cov)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

# grass ~ elevation
ggplot(bada, aes(x = elev_m, y = grass_cov)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()


# is linear model right for moss?
# moss ~ elevation
ggplot(bada, aes(x = elev_m, y = moss_cov)) +
  geom_point() +
  geom_smooth() +
  theme_classic()


# data distributions
hist(bada$heather_cov)
hist(bada$moss_cov)
hist(bada$grass_cov)
hist(bada$lichen_cov)

# transform heather
bada$heather_log <- log(bada$heather_cov)
bada$heather_sqrt <- sqrt(bada$heather_cov)
hist(bada$heather_log)
hist(bada$heather_sqrt)  # probably better
# heather lm
lm_heather <- lm(heather_sqrt ~ elev_m, data = bada)
lm_heather_log <- lm(heather_log ~ elev_m, data = bada)
par(mfrow = c(2,2))
plot(lm_heather)
plot(lm_heather_log)

# transform lichen
bada$lichen_log <- log(bada$lichen_cov)
bada$lichen_sqrt <- sqrt(bada$lichen_cov)
hist(bada$lichen_log)
hist(bada$lichen_sqrt)




bada$moss_log <- log(bada$moss_cov) + 1
hist(bada$moss_log)

chisq.test(bada$elev_m, bada$moss)


## Question 3 - correlation of the two -----------------------------------

# a)
# inv spp richness as a response
hist(meall$inv_order_rich)
hist(log(meall$inv_order_rich))  # looks more normal, therefore I'll use log
meall$inv_rich_log <- log(meall$inv_order_rich)

# log-transformed inv richness
lm_richlog_null <- lm(inv_rich_log ~ 1, data = meall)
lm_richlog_el <- lm(inv_rich_log ~ elev_m, data = meall)
lm_richlog_pl <- lm(inv_rich_log ~ plant_sp_rich, data = meall)
lm_richlog_elpl <- lm(inv_rich_log ~ elev_m + plant_sp_rich, data = meall)
lm_richlog_elplint <- lm(inv_rich_log ~ elev_m * plant_sp_rich, data = meall)
AIC(lm_richlog_null, lm_richlog_el, lm_richlog_pl, lm_richlog_elpl, lm_richlog_elplint)
AICc(lm_richlog_null, lm_richlog_el, lm_richlog_pl, lm_richlog_elpl, lm_richlog_elplint)
# no model was better than the null model. Inv richness doesn't seem to be explained
# by elevation or plant richness
lm_richlog_heather <- lm(inv_rich_log ~ heather_cov, data = meall)
lm_richlog_moss <- lm(inv_rich_log ~ moss_cov, data = meall)
lm_richlog_heatmoss <- lm(inv_rich_log ~ heather_cov + moss_cov, data = meall)
lm_richlog_heatmossint <- lm(inv_rich_log ~ heather_cov * moss_cov, data = meall)
AIC(lm_richlog_null, lm_richlog_heather, lm_richlog_moss, lm_richlog_heatmoss, lm_richlog_heatmossint)
AICc(lm_richlog_null, lm_richlog_heather, lm_richlog_moss, lm_richlog_heatmoss, lm_richlog_heatmossint) # could not find fuction (don't know why)
# vegetation cover seems to explain invertebrate diversity better
# all of the models above are better than the null model
# the best model is for moss coverage as an independent variable
summary(lm_richlog_moss)
plot(inv_rich_log ~ moss_cov, data = meall)
abline(lm_richlog_moss)
par(mfrow = c(2,2))
plot(lm_richlog_moss)  # residuals are not well distributed
# what about not transformed y
lm_rich_moss <- lm(inv_order_rich ~ moss_cov, data = meall)
plot(lm_rich_moss)  # not really better

# test of normality
moss_resid <- resid(lm_richlog_moss)
shapiro.test(moss_resid)  # data are normally distributed

#///////what about poisson instead of log-transformation
glm_rich_null <- glm(inv_order_rich ~ 1, data = meall, family = poisson)
glm_rich_el <- glm(inv_order_rich ~ elev_m, data = meall, family = poisson)
AIC(glm_rich_null, glm_rich_el) # no, the same results as log

glm_rich_heather <- glm(inv_order_rich ~ heather_cov, data = meall, family = poisson)
glm_rich_moss <- glm(inv_order_rich ~ moss_cov, data = meall, family = poisson)
glm_rich_heatmoss <- glm(inv_order_rich ~ heather_cov + moss_cov, data = meall, family = poisson)
glm_rich_heatmossint <- glm(inv_order_rich ~ heather_cov * moss_cov, data = meall, family = poisson)
AIC(glm_rich_null, glm_rich_heather, glm_rich_moss, glm_rich_heatmoss, glm_rich_heatmossint)  # surprisingly, moss is not a good predictor anymore (slightly better)
plot(glm_rich_moss)
summary(glm_rich_moss)

#///////log transformation or poisson? log-transformation gives more interesting results.

# conclusions: Species richness of invertebrates seems to be affected by percentage 
# cover of moss (p = 0.013) if we log-transform richness. Residuals however are not evenly distributed.
# If we use poisson distribution instead, no model is better than the null model.
# Either way, inv richness was not affected by elevation alone.
par(mfrow = c(1,1))



# b)
# inv diversity as a response
hist(meall$shannon_diversity) # looks ok

lm_div_null <- lm(shannon_diversity ~ 1, data = meall)
lm_div_el <- lm(shannon_diversity ~ elev_m, data = meall)
lm_div_pl <- lm(shannon_diversity ~ plant_sp_rich, data = meall)
lm_div_elpl <- lm(shannon_diversity ~ elev_m + plant_sp_rich, data = meall)
lm_div_elplint <- lm(shannon_diversity ~ elev_m * plant_sp_rich, data = meall)
AIC(lm_div_null, lm_div_el, lm_div_pl, lm_div_elpl, lm_div_elplint)
AICc(lm_div_null, lm_div_el, lm_div_pl, lm_div_elpl, lm_div_elplint)
# no model is better than the null model

lm_div_heather <- lm(shannon_diversity ~ heather_cov, data = meall)
lm_div_moss <- lm(shannon_diversity ~ moss_cov, data = meall)
lm_div_heatmoss <- lm(shannon_diversity ~ heather_cov + moss_cov, data = meall)
lm_div_heatmossint <- lm(shannon_diversity ~ heather_cov * moss_cov, data = meall)
AIC(lm_div_null, lm_div_heather, lm_div_moss, lm_div_heatmoss, lm_div_heatmossint)
# Inv diversity seems to be best explained by heather and moss with interaction
# Second best model is explained by moss only.
# If AICc would work, maybe moss would be the best model.

# moss and heather with interaction
summary(lm_div_heatmossint)
par(mfrow = c(2,2))
plot(lm_div_heatmossint)  # one outlier - plot 1

# moss only
summary(lm_div_moss)
plot(lm_div_moss)  # not horrible, better than the one above, so I would choose this one
par(mfrow = c(1,1))  
plot(shannon_diversity ~ moss_cov, data = meall)
abline(lm_div_moss)

# conclusions: moss coverage seems to be the best explanation for inv diversity

# c)
# inv abundance as a response
hist(meall$inv_abun) # right-skewed, and they are also counts, so would go for glm poisson

glm_abu_null <- glm(inv_abun ~ 1, data = meall, family = poisson)
glm_abu_el <- glm(inv_abun ~ elev_m, data = meall, family = poisson)
glm_abu_pl <- glm(inv_abun ~ plant_sp_rich, data = meall, family = poisson)
glm_abu_elpl <- glm(inv_abun ~ elev_m + plant_sp_rich, data = meall, family = poisson)
glm_abu_elplint <- glm(inv_abun ~ elev_m * plant_sp_rich, data = meall, family = poisson)
AIC(glm_abu_null, glm_abu_el, glm_abu_pl, glm_abu_elpl, glm_abu_elplint)
AICc(glm_abu_null, glm_abu_el, glm_abu_pl, glm_abu_elpl, glm_abu_elplint)
# abundance is best explained by elevation and plant richness without interaction
# but elevation only is close 2nd
# let's check the two models
par(mfrow = c(2,2))  
plot(glm_abu_elpl)  # this has two strong outliers, better normality
plot(glm_abu_el)  # this has one strong outlier, better residuals

par(mfrow = c(1,1))  
summary(glm_abu_elpl)
summary(glm_abu_el)
plot(inv_abun ~ elev_m, data = meall)
abline(glm_abu_el)  # strange trendline

glm_abu_heather <- glm(inv_abun ~ heather_cov, data = meall, family = poisson)
glm_abu_moss <- glm(inv_abun ~ moss_cov, data = meall, family = poisson)
glm_abu_heatmoss <- glm(inv_abun ~ heather_cov + moss_cov, data = meall, family = poisson)
glm_abu_heatmossint <- glm(inv_abun ~ heather_cov * moss_cov, data = meall, family = poisson)
glm_abu_elmoss <- glm(inv_abun ~ elev_m + moss_cov, data = meall, family = poisson)
glm_abu_elmossint <- glm(inv_abun ~ elev_m * moss_cov, data = meall, family = poisson)
AIC(glm_abu_null, glm_abu_heather, glm_abu_moss, glm_abu_heatmoss, glm_abu_heatmossint,
    glm_abu_elmoss, glm_abu_elmossint, glm_abu_el)
# all are better than null. Best is moss and elevation with interaction.

par(mfrow = c(2,2)) 
plot(glm_abu_moss)  # residuals not perfect, neither normality, but not too bad
plot(glm_abu_elmossint)  # looks worse than just moss, so I would go for moss
par(mfrow = c(1,1)) 
plot(inv_abun ~ moss_cov, data = meall)
abline(glm_abu_moss)  # not a bad graph, higher variability in high moss coverage

summary(glm_abu_moss)

## conlcusions: abundance is quite well explained by elevation and moss cover (separately).
# however, moss cover is a better predictor. 



## NMDS --------------------------------------------------------------

library(vegan)
library(tidyverse)

inv <- meall[,16:26]
# make NMDS
inv %>%
  metaMDS(trace = F) %>%
  ordiplot(type = "none") %>%
  text("sites")

# First step is to calculate a distance matrix. See PCOA for more information about the distance measures
# Here we use bray-curtis distance, which is recommended for abundance data
dist <- vegdist(inv,  method = "bray")

# In this part, we define a function NMDS.scree() that automatically 
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Use the function that we just defined to choose the optimal nr of dimensions
NMDS.scree(dist)

# Because the final result depends on the initial 
# random placement of the points 
# we`ll set a seed to make the results reproducible
set.seed(2)

# Here, we perform the final analysis and check the result
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F)
# Do you know what the trymax = 100 and trace = F means?
# Let's check the results
NMDS1

# If you don`t provide a dissimilarity matrix, metaMDS automatically applies Bray-Curtis. So in our case, the results would have to be the same
NMDS2 <- metaMDS(inv, k = 2, trymax = 100, trace = F, autotransform = FALSE)
NMDS2

# Let’s check the results of NMDS1 with a stressplot
stressplot(NMDS1)

plot(NMDS1, type = "t")
# There are no species scores (same problem as we encountered with PCoA).
# We can work around this problem, by giving metaMDS the original community matrix
# as input and specifying the distance measure.
NMDS3 <- metaMDS(inv, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
plot(NMDS3)
plot(NMDS3, display = "sites", type = "n")
points(NMDS3, display = "sites", col = "red", cex = 1.25)
text(NMDS3, display ="species")

# Alternatively, you can use the functions ordiplot and orditorp
ordiplot(NMDS3, type = "n")
orditorp(NMDS3, display = "species", col = "red", air = 0.01)
orditorp(NMDS3, display = "sites", cex = 1.1, air = 0.01)

# We now have a nice ordination plot and we know which plots have a similar species composition.
# The next question is: Which environmental variable is driving the observed differences
# in species composition? We can do that by correlating environmental variables with our ordination axes.
# Therefore, we will use a second dataset with environmental variables (sample by environmental variables).
# We continue using the results of the NMDS.

# load the second dataset
env <- meall %>% 
  select(elev_m, plant_sp_rich, heather_cov, moss_cov, grass_cov, lichen_cov)
head(env)
names(env) <- c("elevation", "plant richness", "heather %", "moss %", "grass %", "lichen %")

# The function envfit will add the environmental variables as vectors to the ordination plot
ef <- envfit(NMDS3, env, permu = 999)
ef

# The two last columns are of interest: the squared correlation coefficient and the associated p-value
# Plot the vectors of the significant correlations and interpret the plot
plot(NMDS3, type = "t", display = "sites")
plot(ef, p.max = 0.80)  # should be p.max = 0.05
