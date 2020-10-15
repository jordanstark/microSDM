# NicheMapR for Landscape Ecology SDM project

#### setup ####
#not a CRAN package, installed from github with devtools::install_github('mrke/NicheMapR')
 # you may need to run R as an admin for that!
 # I also needed to install package 'futile logger'
 # also install global monthly climate means (this is used for the solar radiation calculation if nothing else))
 #get.global.climate(folder="C:/Users/Jordan/Desktop/Landscape ecology/Term_project/Global_Climate")
 # also requires futile.logger package which is on CRAN so can be installed normally

# for the finest scale models, this integrates with microclima
#devtools::install_github("ilyamaclean/microclima")

# directories
sp_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/Species_data/"

#packages
library(NicheMapR)
library(microclima)
library(RNCEP) #to download climate data 
library(lubridate) # deals with dates
library(ggplot2)
library(readr)
library(sp)
library(raster)

# import & clean species data
pls <- read_delim(paste(sp_path,"Pinkladyslip.csv",sep=""),delim="\t",col_names = TRUE, na=c("","NA"))
pls_simple <- pls[,c("scientificName","stateProvince","occurrenceStatus","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","eventDate")]
pls_good <- pls_simple[pls_simple$coordinateUncertaintyInMeters<100 & !is.na(pls_simple$coordinateUncertaintyInMeters),] #2878 observations
pls_ny <- pls_good[pls_good$stateProvince=="New York",] # 390 observations

coords <- pls_ny[,c("decimalLongitude","decimalLatitude")] #could substitute in pls_good for all pts

coords <- coords[complete.cases(coords),]

pls_spatial <- SpatialPoints(coords,proj4string=CRS("+proj=longlat"))

plot(pls_spatial)


#### better version using microclima ####

# version for 30m2 microclimate (I think this is what we want to use eventually)

r <- raster(pls_spatial)


dem <- get_dem(r=r,resolution=30)
  # this automatically downloads a DEM that can be used for all microclimate models
  # could be saved using write.raster() - not sure if this would be faster than the auto download below
  # 'dims' are number of cells
  # resolution in m

save = 2 # save=1 to save input data; save=2 to use saved data

new_micro <- micro_ncep(loc=coords[1], dstart=dstart, dfinish=dfinish,ERR=3,save=save)
#Err is tolerance for soil temp; at default 1.5 the model is crashing
# this runs in ~ 13 seconds if using save=2
# It is possible to do this for grids, 
# but if I am understanding the SDMs right we actually just need it for presence and pseudo-absence points

new_temp <- new_micro$microclima.out$hourlydata
new_temp$date <- date(new_temp$obs_time)

daily_temp <- aggregate(temperature ~ date, new_temp, mean)
daily_precip <- new_micro$microclima.out$dailyprecip
  # I believe this is just interpolated but may be better than nothing

new_moist <- as.data.frame(new_micro$soilmoist)
  # not sure if this is actually affected by using microclima package
new_moist$date <- date(new_micro$dates)
daily_surfacemoist <- aggregate(WC0cm ~ date, new_moist, mean) 


# compare outputs

soilmoist$date <- date(micro$dates)
usasoilmoist <- aggregate(WC0cm ~ date, soilmoist, mean) 

metout$date <- date(micro$dates)
usatemp <- aggregate(TAREF ~ date, metout, mean)

all_dat <- data.frame(usa_temp = usatemp$TAREF,
                      usa_precip = micro$RAINFALL,
                      usa_surfacemoist = usasoilmoist$WC0cm,
                      new_temp = daily_temp$temperature,
                      new_precip = daily_precip,
                      new_surfacemoist = daily_surfacemoist$WC0cm,
                      date = daily_temp$date)

ggplot(all_dat, aes(x=date)) +
  geom_line(aes(y=usa_temp),color="blue") +
  geom_line(aes(y=new_temp),color="green") +
  theme_classic()

ggplot(all_dat, aes(x=date)) +
  geom_line(aes(y=usa_precip),color="blue") +
  geom_line(aes(y=new_precip),color="green") +
  theme_classic()

ggplot(all_dat, aes(x=date)) +
  geom_line(aes(y=usa_surfacemoist),color="blue") +
  geom_line(aes(y=new_surfacemoist),color="green") +
  theme_classic()

# so, they are very different (especially moisture) -- but it's possible this is just due to different assumptions


# more detailed treatment & all parameters described here: https://mrke.github.io/NicheMapR/inst/doc/microclimate-IO
# also in the supplement to Kearney et al 2020 which describes micro_ncep

