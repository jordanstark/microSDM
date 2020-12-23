# figures and analysis of microclimate and SDMS
# landscape ecology, fall 2020


#### setup ####

library(ggplot2)
library(lubridate)
library(sp)
library(raster)
library(tmap)
library(tmaptools)

sp_path <- "E:/LandscapeEco/Species_data/"
gis_path <- "E:/LandscapeEco/GIS/"
climate_path <- "E:/LandscapeEco/NicheMapR/"

states <- getData("GADM",country="USA",level=1) #download state data
NY <- states[states$NAME_1=="New York",]
NY_poly <- SpatialPolygons(NY@polygons,proj4string=crs("+proj=longlat"))

#### climate across the state figure ####

bioclim <- read.csv(paste(sp_path,"bioclim_20ymean_shademicroclimate.csv",sep=""))
bioclim$X <- NA

bioclim_sp <- SpatialPointsDataFrame(coords=cbind(bioclim$lon,bioclim$lat),proj4string=crs("+proj=longlat"),data=bioclim)

ggplot(bioclim,aes(x=lon,y=lat,color=bio1))+
  geom_point()

tm_shape(NY_poly) +
  tm_polygons(alpha=0.2) +
tm_shape(bioclim_sp) +
  tm_dots(col="bio1",size=1,palette="-RdBu",shape="pls_presence",shapes=c(15,19),
          legend.show=F,title.shape="Pink Lady Slipper Presence")

tm_shape(NY_poly) +
  tm_polygons(alpha=0.2) +
tm_shape(bioclim_sp) +
  tm_dots(col="bio2",size=1,palette="Greys",shape="pls_presence",shapes=c(15,19),
          legend.show=F,title.shape="Pink Lady Slipper Presence")

#### climate differences between presence/absence figure ####

library(tidyr)

bioclim_long <- pivot_longer(bioclim,cols=starts_with("bio"),values_to="value",names_to="variable")

ggplot(bioclim_long,aes(x=pls_presence,y=value)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free") +
  theme_classic()

#### compare with and without shading ####
no_shade <- read.csv(paste(sp_path,"bioclim_20ymean_simplemicroclimate.csv",sep=""))
no_shade$X <- NULL
no_shade_long <- pivot_longer(no_shade,cols=starts_with("bio"),values_to="value_noshade",names_to="variable")

all <- merge(bioclim_long, no_shade_long)
all$shade_diff <- all$value - all$value_noshade

ggplot(all,aes(x=pls_presence,y=shade_diff)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free") +
  theme_classic()

ggplot(all[all$variable!="bio4",],aes(x=variable,y=shade_diff)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") + ylab("shade climate - unshaded climate")


#### timeseries ####
raw_climate_dat <- read.csv(paste(climate_path,"20y_Shade_microclimate_summary.csv",sep=""))
raw_climate_dat$X <- NULL
raw_climate_dat$date <- ymd(raw_climate_dat$date)

ggplot(raw_climate_dat, aes(x=date,y=maxt,color=pls_presence)) +
  geom_line(aes(group=ptID),alpha=0.1) +
  theme_classic() +
  ylab("") + xlab("") 

#### compare bioclim and nichemapr climates ####
bioclim_1km <- getData('worldclim',var="bio",res=0.5,lat=42,lon=-75,)

