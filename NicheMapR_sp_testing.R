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
gis_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/GIS/elevation/"

#packages
library(NicheMapR)
library(microclima)
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

pls_spatial <- SpatialPointsDataFrame(coords,proj4string=CRS("+proj=longlat"),data=coords)

plot(pls_spatial)


dem1 <- raster(paste(gis_path,"ned30m40073.tif",sep=""))

pls_utm <- spTransform(pls_spatial,crs(dem1))

pls_dem1 <- crop(pls_utm,dem1)
dem1coords <- pls_dem1@data

# calculating inputs for NicheMapR with 30m DEM
slope <- terrain(dem1,opt="slope",unit="degrees")
aspect <- terrain(dem1,opt="aspect",unit="degrees")

# extract for points
pls_dem1$elev <- extract(dem1,pls_dem1)
pls_dem1$slope <- extract(slope,pls_dem1)
pls_dem1$aspect <- extract(aspect,pls_dem1)

dem1data <- data.frame(pls_dem1@data[,c("elev","slope","aspect")])

# model parameters
runshade <- 0 #run for only one shade value
runmoist <- 1 #run soil moisture model
soilgrids <- 1 #get soil properties from soilgrids.org
                  # this increases time from ~18 sec/site to ~32 sec/site



#### additional data/params that would be good to include
#pre-downloading met data; ideally would use same data as macro climate models
#Usrhyt; heigh for temp, humidity estimation - default is 0.01m, consider using different for different orgs?
#minshade; % shade used for calcs - try to get an estimate of canopy cover? This may also need to vary seasonally
# hori - horizon angles; can calculate with horizon package but it takes a while and probably not accurate near edge of DEM
# LAI - leaf area index, for calculating transpiration effect on soil moisture
# snow model !!!




# function to extract data of interest from model result
extractdailydat <- function(micro_out) {
  soildf <- data.frame(micro_out$soilmoist)
  metdf <- data.frame(micro_out$metout)
  dates <- micro_out$dates
  
  alldf_6h <- cbind(soildf,metdf,dates)
  alldf_6h$date <- date(alldf_6h$dates)
  
  daily_summary <- aggregate(cbind(WC0cm,WC10cm,WC100cm,TALOC,RHLOC) ~ date, alldf_6h, mean)
  names(daily_summary) <- c("date","moist_0","moist_10","moist_100","AirT","Rel.Hum")
  daily_summary$Precip <- micro_out$RAINFALL
  
  return(daily_summary)
  
} 

summarylist <- list()

for(i in 1:length(dem1coords$decimalLongitude)) {
  ptdat <- micro_usa(loc=as.numeric(dem1coords[i,]), dstart="01/01/2019", dfinish="31/12/2019",
                     slope=dem1data$slope[i],elev=dem1data$elev[i],aspect=dem1data$aspect[i],
                     runshade=runshade,runmoist=runmoist,soilgrids=soilgrids)
  summarylist[[i]] <- extractdailydat(ptdat)
  summarylist[[i]]$pt <- i
}

summarydf <- do.call(rbind,summarylist)
summarydf$pt <- factor(summarydf$pt)



ggplot(summarydf, aes(x=date,y=AirT,color=pt)) +
  geom_line() +
  theme_classic()
ggplot(summarydf, aes(x=date,y=moist_0,color=pt)) +
  geom_line() +
  theme_classic()
ggplot(summarydf, aes(x=date,y=Rel.Hum,color=pt)) +
  geom_line() +
  theme_classic()
ggplot(summarydf, aes(x=date,y=Precip,color=pt)) +
  geom_line() +
  theme_classic()





