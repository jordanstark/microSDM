# Prep for NicheMapR models for NYS points
# Landscape ecology final project
# Jordan Stark

#### setup ####
# paths
gis_path <- "E:/LandscapeEco/GIS/"

# packages
library(readr) # to read species data
library(sp)
library(raster)
library(stringr)
#library(horizon) # to calculate sky view or horizon angle (not currently working)


#### generate points ####
  # select points in NY for background data
  states <- getData("GADM",country="USA",level=1) #download state data
  NY <- states[states$NAME_1=="New York",]
  NY_poly <- SpatialPolygons(NY@polygons)
  crs(NY_poly) <- crs(NY)
  
  set.seed(9842368)
  
  background <- spsample(NY_poly,10000,type="regular") #10,000 points is way larger than 30m but starting somewhere...
  background_spdf <- SpatialPointsDataFrame(background,data=data.frame(background@coords))
  names(background_spdf@data) <- c("decimalLongitude","decimalLatitude")
  crs(background_spdf) <- proj4string(CRS("+proj=longlat"))
  
  write.csv(background_spdf,paste(gis_path,"NY_points.csv",sep=""))

  
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
  
  
# prep fPAR files for shading
  # fpar_files <-  list.files(paste(gis_path,"FPAR/FPAR/",sep=""),full.names=T)
  # fpar_stack <- stack(fpar_files)
  # 
  # qual_files <- list.files(paste(gis_path,"FPAR/QC/",sep=""),full.names=T)
  # qual_stack <- stack(qual_files)
  # 
  # #qc_lookup <- read.csv(paste(gis_path,"FPAR/metadata/MCD15A3H-006-FparLai-QC-lookup.csv",sep=""))
  #   # probably could automate use of this but I'm manually choosing good #s for now
  # 
  # good_qc_vals <- c(0,2,16,18,32,34,48,50) # not significant clouds, good quality, main algorithm
  # 
  # 
  # for(i in 1:nlayers(fpar_stack)){
  #   dat <- fpar_stack[[i]]
  #   doy <- str_split(names(fpar_stack[[i]]),pattern="_",simplify=T)[,4]
  #   
  #   qc <- qual_stack[[i]]
  #   qc[qc %in% good_qc_vals] <- 1
  #   good_dat <- mask(dat,qc,inverse=T,maskvalue=1,filename=paste(gis_path,"FPAR/CleanFPAR/",doy,sep=""))
  # } #start at9:24 on 11/3, ~29days/min with 1562 layers so should take an hour
  # 
  # rm(qual_stack)
  # rm(fpar_stack) # these are quite large so R runs better without them
  
  clean_fpar_files <- list.files(paste(gis_path,"FPAR/CleanFPAR/",sep=""),pattern=".gri",full.names=T)
  clean_fpar <- stack(clean_fpar_files)
  
  background_pts <- read.csv(paste(gis_path,"NY_points.csv",sep=""))
  background_spdf <- SpatialPointsDataFrame(coords=background_pts[,c("decimalLongitude","decimalLatitude")],
                                            proj4string=crs("+proj=longlat"),data=background_pts[,c("decimalLongitude","decimalLatitude")])
  background_pts <- NULL
  
  background_spdf$ID <- 1:length(background_spdf$decimalLongitude)
  pts_sinu <- spTransform(background_spdf,crs(clean_fpar))
  
  
  fpar_points <- data.frame(raster::extract(clean_fpar,pts_sinu))
  pointsOnly <- fpar_points # for rowSums
  fpar_points$decimalLatitude  <- background_spdf$decimalLatitude
  fpar_points$decimalLongitude <- background_spdf$decimalLongitude
  fpar_points$ID <- background_spdf$ID
  ID_cols <- c("decimalLatitude","decimalLongitude","ID")
  

  
  fpar_points_clean <- fpar_points[rowSums(pointsOnly,na.rm=T)>0,] #remove 1271 points with no obs
  
  library(tidyr)
  library(lubridate)

  fpar_long <- pivot_longer(fpar_points_clean,cols=starts_with("MCD"),
                            names_to="filename",values_to="FPAR")
  fpar_long$doy <- str_split(fpar_long$filename,pattern="_",simplify=T)[,4]
  fpar_long$date <- str_split(fpar_long$doy,pattern="doy",simplify=T)[,2]
  fpar_long$date <- parse_date_time(fpar_long$date,orders="Yj")

  
  # add all missing days without records
  all_dates <- seq(min(fpar_long$date),max(fpar_long$date),by="days")
  all_combos <- data.frame(expand.grid(all_dates,unique(fpar_long$ID)))
  names(all_combos) <- c("date","ID")

  
  fpar_long <- fpar_long[,c("decimalLatitude","decimalLongitude","ID","FPAR","date")]
  fpar_long <- data.frame(fpar_long)
  
  rm(clean_fpar,fpar_points,NY,NY_poly,states,pts_sinu,background_spdf,background)

  
  fpar_pad <- merge(fpar_long,all_combos,by=c("date","ID"),all.y=T)

  all_combos[,c("ID")] <- NA
  
  fpar_pad$year <- year(fpar_pad$date)
  fpar_pad$day <- yday(fpar_pad$date)  
  
  fpar_pad <- fpar_pad[order(fpar_pad$date),] 
  
  fpar_list <- split(fpar_pad,fpar_pad$ID)
  # fill missing values with previous value (might want to consider something else but this is simple)
  avg_fpar_list <- list()
  
  for(i in 1:length(fpar_list)) {
    #fpar_list[[i]]$decimalLatitude <- fpar_list[[i]]$decimalLatitude[1]
    #fpar_list[[i]]$decimalLongitude <- fpar_list[[i]]$decimalLongitude[1]

    dat <- fill(fpar_list[[i]],FPAR,.direction="downup") #fills down except the first few rows which are filled up
    avg_fpar_list[[i]] <- aggregate(FPAR ~ day,dat,mean)
    
    avg_fpar_list[[i]]$decimalLatitude <- fpar_list[[i]]$decimalLatitude[1]
    avg_fpar_list[[i]]$decimalLongitude <- fpar_list[[i]]$decimalLongitude[1]
    avg_fpar_list[[i]]$ID <- fpar_list[[i]]$ID[1]
    
    
   }
  
  avg_fpar <- do.call(rbind,avg_fpar_list)

  
  # plot to check
  library(ggplot2)
  ggplot(avg_fpar,aes(x=day,y=FPAR,color=factor(ID))) +
    geom_line() +
    theme_classic() +
    theme(legend.position="none")
  
  write.csv(avg_fpar,paste(gis_path,"NY_FPAR.csv",sep=""),row.names=F)
  
  
  # save new point list including only these points
  good_pts <- unique(avg_fpar[,c("decimalLatitude","decimalLongitude")])
  write.csv(good_pts,paste(gis_path,"NY_points_good.csv",sep=""),row.names=F)


  
  