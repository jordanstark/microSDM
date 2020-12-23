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
#sp_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/Species_data/"
#gis_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/GIS/elevation/"
sp_path <- "E:/LandscapeEco/Species_data/"
gis_path <- "E:/LandscapeEco/GIS/elevation"
out_path <- "E:/LandscapeEco/NicheMapR/"


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

# locate DEMs
demlist <- list.files(gis_path,pattern="tif$",full.names=T)



#### for each DEM, import, crop points, and calculate params for microclimate ####


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

start_time <- Sys.time()
all_pt_list <- list()

for(i in 1:length(demlist)) {
  #import DEM and extract associated points
  dem <- raster(demlist[i])
  pls_utm <- spTransform(pls_spatial, crs(dem)) #transform points to crs of DEM
  
  
  if(length(crop(pls_utm,dem)) != 0){
    pts_in_dem <- crop(pls_utm,dem)
    check_in_dem <- extract(dem,pts_in_dem)
    pts_in_dem <- pts_in_dem[!is.na(check_in_dem),] 
      #this removes pts that are in the bounding box but not actually on the DEM
    
    new_ext <- extent(pts_in_dem) + 1000
    dem <- crop(dem,new_ext)
      # this speeds up calculation by only including area within 1km of max point spread
    
    
 
  
    coords <- pts_in_dem@data
    
    #calculate other needed values
    slope <- terrain(dem,opt="slope",unit="degrees")
    aspect <- terrain(dem,opt="aspect",unit="degrees")
    
    #extract values at each point
    pts_in_dem$elev <- extract(dem,pts_in_dem)
    pts_in_dem$slope <- extract(slope,pts_in_dem)
    pts_in_dem$aspect <- extract(aspect,pts_in_dem)
    
    ptdata <- data.frame(pts_in_dem@data[,c("elev","slope","aspect")])
    
    #set model parameters
    runshade <- 0 #run for only one shade value
    runmoist <- 1 #run soil moisture model
    soilgrids <- 1 #get soil properties from soilgrids.org
        # this increases time from ~18 sec/site to ~32 sec/site
    ERR <- 2 # default 1.5 leads to some sites failing
                        
    
    dem_pt_list <- list()
    
    for(j in 1:length(coords$decimalLongitude)) {
      skip <- FALSE
      
      dat <- tryCatch(micro_usa(loc=as.numeric(coords[j,]), dstart="01/01/2019", dfinish="31/12/2019",
                         slope=ptdata$slope[j],elev=ptdata$elev[j],aspect=ptdata$aspect[j],
                         runshade=runshade,runmoist=runmoist,soilgrids=soilgrids),
                      error=function(e){skip <<- TRUE})
      
      if(skip==F) {
        dem_pt_list[[j]] <- extractdailydat(dat) 
        
        write.csv(dat,paste(out_path,"microusa_2019_",i,"_",j,".csv",sep=""))
        
      } else {
        dem_pt_list[[j]] <- data.frame(date=NA, moist_0=NA, moist_10=NA, moist_100=NA,
                                       AirT=NA, Rel.Hum=NA, Precip=NA)
        }

      dem_pt_list[[j]]$lat <- as.numeric(coords[j,2])
      dem_pt_list[[j]]$lon <- as.numeric(coords[j,1])
      dem_pt_list[[j]]$pt_dem <- j 
    }
    
    dem_pt_df <- do.call(rbind,dem_pt_list)
    
  } else { dem_pt_df <- data.frame(date=NA, moist_0=NA, moist_10=NA, moist_100=NA,
                                   AirT=NA, Rel.Hum=NA, Precip=NA,
                                   lat=NA, lon=NA, pt_dem=NA)}
  
  all_pt_list[[i]] <- dem_pt_df
  all_pt_list[[i]]$dem_num <- i
  
}


all_pt_dat <- do.call(rbind,all_pt_list)
all_pt_dat$ptID <- paste(all_pt_dat$pt,all_pt_dat$pt_dem,sep=".")

end_time <- Sys.time()

start_time
end_time




ggplot(all_pt_dat, aes(x=date,y=AirT,color=ptID)) +
  geom_line(alpha=0.2) +
  theme_classic() +
  theme(legend.position="none")
ggplot(all_pt_dat, aes(x=date,y=moist_0,color=ptID)) +
  geom_line(alpha=0.2) +
  theme_classic() +
  theme(legend.position="none")
ggplot(all_pt_dat, aes(x=date,y=Rel.Hum,color=ptID)) +
  geom_line(alpha=0.2) +
  theme_classic() +
  theme(legend.position="none")
ggplot(all_pt_dat, aes(x=date,y=Precip,color=ptID)) +
  geom_line(alpha=0.2) +
  theme_classic() +
  theme(legend.position="none")



annual_mean <- aggregate(cbind(moist_0,AirT,Precip) ~ ptID + lat + lon, 
                         all_pt_dat, mean)
annual_min <- aggregate(cbind(moist_0,AirT,Precip) ~ ptID + lat + lon, 
                         all_pt_dat, min)
annual_max <- aggregate(cbind(moist_0,AirT,Precip) ~ ptID + lat + lon, 
                         all_pt_dat, max)
annual_sd <- aggregate(cbind(moist_0,AirT,Precip) ~ ptID + lat + lon, 
                        all_pt_dat, sd)

mean_SPDF <- SpatialPointsDataFrame(coords=annual_mean[,c("lon","lat")],
                                    data=annual_mean,
                                    proj4string=CRS("+proj=longlat"))

min_SPDF <- SpatialPointsDataFrame(coords=annual_min[,c("lon","lat")],
                                    data=annual_min,
                                    proj4string=CRS("+proj=longlat"))

max_SPDF <- SpatialPointsDataFrame(coords=annual_max[,c("lon","lat")],
                                    data=annual_max,
                                    proj4string=CRS("+proj=longlat"))

sd_SPDF <- SpatialPointsDataFrame(coords=annual_sd[,c("lon","lat")],
                                    data=annual_sd,
                                    proj4string=CRS("+proj=longlat"))


plot(sd_SPDF,col=mean_SPDF$AirT)

ggplot(annual_min, aes(x=lon, y=lat, color=AirT)) +
  geom_point()
