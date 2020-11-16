# Prep for NicheMapR models for pink lady slipper points
# Landscape ecology final project
# Jordan Stark

# inputs
  # species data from GBIF 
    # "Pinkladyslip.csv", in folder sp_path
  # state boundary data - 
    #automatically downloaded the first time this script runs
  # 30m DEM including all of NYS 
    # "Full_DEM/NASADEM_NC.001_NASADEM_HGT_doy2000042_aid0001.tif"
    # from AppEEARS NASADEM, in folder gis_path

# outputs
  # cleaned points with pink lady slipper and background points 
    # "all_pls_points.csv", in folder sp_path
  # raster of slope in degrees
    # "NYS_slope_degrees.gri" in folder gis_path
  # raster of aspect in degrees
    # "NYS_aspect_degrees.gri" in folder gis_path
  

#### setup ####
# paths
sp_path <- "E:/LandscapeEco/Species_data/"
gis_path <- "E:/LandscapeEco/GIS/"

# packages
library(readr) # to read species data
library(sp)
library(raster)
#library(horizon) # to calculate sky view or horizon angle (not currently working)


#### clean species data and generate background points ####

# import & clean species data
  pls <- read_delim(paste(sp_path,"Pinkladyslip.csv",sep=""),delim="\t",col_names = TRUE, na=c("","NA"))
  pls_simple <- pls[,c("scientificName","stateProvince","occurrenceStatus","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","eventDate")]
  pls_good <- pls_simple[pls_simple$coordinateUncertaintyInMeters<100 & !is.na(pls_simple$coordinateUncertaintyInMeters),] #2878 observations
  pls_ny <- pls_good[pls_good$stateProvince=="New York",] # 390 observations
                                                          # with 30m accuracy, 191 obs ... so may be possible?
  
  coords <- pls_ny[,c("decimalLongitude","decimalLatitude")] #could substitute in pls_good for all pts
  
  coords <- coords[complete.cases(coords),]



# remove points that are too close together (based on https://www.molecularecologist.com/2013/04/species-distribution-models-in-r/)
  longrid <- seq(min(coords$decimalLongitude),max(coords$decimalLongitude),0.05) # for ~5km grids, could do smaller
  latgrid <- seq(min(coords$decimalLatitude),max(coords$decimalLatitude),0.05)
  
  coords_sub <- data.frame(coords)
  
  subset <- c()
  for(i in 1:(length(longrid)-1)) {
    for(j in 1:(length(latgrid)-1)){
      gridsq <- subset(coords_sub, decimalLatitude > latgrid[j] & decimalLatitude < latgrid[j+1] & decimalLongitude > longrid[i] & decimalLongitude < longrid[i+1])
      if(dim(gridsq)[1]>0){
        subset <- rbind(subset, gridsq[sample(1:dim(gridsq)[1],1 ), ])
  
      }
  
    }
  
  }
  
  pls_points <- SpatialPointsDataFrame(subset,proj4string=CRS("+proj=longlat"),data=subset) # 109 points
  pls_points$pls_presence <- T




# select random points in NY for background data
  states <- getData("GADM",country="USA",level=1) #download state data
  NY <- states[states$NAME_1=="New York",]
  NY_poly <- SpatialPolygons(NY@polygons)
  
  set.seed(68453)
  
  background <- spsample(NY_poly,length(pls_points$decimalLongitude),type="random")
  background_spdf <- SpatialPointsDataFrame(background,data=data.frame(background@coords))
  names(background_spdf@data) <- c("decimalLongitude","decimalLatitude")
  background_spdf$pls_presence <- F
  crs(background_spdf) <- proj4string(CRS("+proj=longlat"))

# combine background and data points

  
  pls_spatial <- rbind(pls_points,background_spdf)
  plot(pls_spatial)
  write.csv(pls_spatial@data,paste(sp_path,"all_pls_points.csv",sep=""),row.names=F)
  
#### Calculate GIS variables needed for microclimate model ####
#import DEM
  DEM <- raster(paste(gis_path,"Full_DEM/NASADEM_NC.001_NASADEM_HGT_doy2000042_aid0001.tif",sep=""))
  
# calculate other variables from DEM - this takes a while
  slope <- terrain(DEM,opt="slope",unit="degrees")
  writeRaster(slope,paste(gis_path,"NYS_slope_degrees",sep=""))

  aspect <- terrain(DEM,opt="aspect",unit="degrees")
  writeRaster(aspect,paste(gis_path,"NYS_aspect_degrees",sep=""))

  #svf <- svf(DEM) #sky view factor
  #writeRaster(svf,paste(gis_path,"NYS_svf",sep="")) # not using at the moment
  
# horizon angle - I have not actually run this!
  
  # #### This takes AT LEAST an hour per raster so I haven't run yet
  # 
  # 
  # start_time <- Sys.time()
  # h_0 <- horizonSearch(DEM,azimuth=0)
  # end_time <- Sys.time()
  # start_time
  # end_time
  # 
  # h_15 <- horizonSearch(DEM,azimuth=15)
  # h_30 <- horizonSearch(DEM,azimuth=30)
  # h_45 <- horizonSearch(DEM,azimuth=45)
  # h_60 <- horizonSearch(DEM,azimuth=60)
  # h_75 <- horizonSearch(DEM,azimuth=75)
  # h_90 <- horizonSearch(DEM,azimuth=90)
  # h_105 <- horizonSearch(DEM,azimuth=105)
  # h_120 <- horizonSearch(DEM,azimuth=120)
  # h_135 <- horizonSearch(DEM,azimuth=135)
  # h_150 <- horizonSearch(DEM,azimuth=150)
  # h_165 <- horizonSearch(DEM,azimuth=165)
  # h_180 <- horizonSearch(DEM,azimuth=180)
  # h_195 <- horizonSearch(DEM,azimuth=195)
  # h_210 <- horizonSearch(DEM,azimuth=210)
  # h_225 <- horizonSearch(DEM,azimuth=225)
  # h_240 <- horizonSearch(DEM,azimuth=240)
  # h_255 <- horizonSearch(DEM,azimuth=255)
  # h_270 <- horizonSearch(DEM,azimuth=270)
  # h_285 <- horizonSearch(DEM,azimuth=285)
  # h_300 <- horizonSearch(DEM,azimuth=300)
  # h_315 <- horizonSearch(DEM,azimuth=315)
  # h_330 <- horizonSearch(DEM,azimuth=330)
  # h_345 <- horizonSearch(DEM,azimuth=345)
  # 
  # h_stack <- stack(h_0,h_15,h_30,h_45,h_60,h_75,
  #                  h_90,h_105,h_120,h_135,h_150,h_165,
  #                  h_180,h_195,h_210,h_225,h_240,h_255,
  #                  h_270,h_285,h_300,h_315,h_330,h_345)
  # 
  # writeRaster(h_stack,paste(gis_path,"horizon_angle_stack",sep=""))
  # 