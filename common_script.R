# Project 7 - How is invertebrate composition affected by elevation in Glenmore?
# Field course - September 2023
# Common script
# Adela, Adela, Emily, Lubomir, Nina

library(tidyverse)

bada <- read.csv("meall_data.csv")


## Exploration of vegetation ---------------------------------------------

# plant richness ~ elevation
ggplot(bada, aes(x = elev_m, y = plant_sp_rich)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

# heather ~ elevation
ggplot(bada, aes(x = elev_m, y = heather_cov)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

# moss ~ elevation
ggplot(bada, aes(x = elev_m, y = moss_cov)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

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