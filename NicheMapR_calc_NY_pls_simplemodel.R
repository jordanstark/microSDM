# NicheMapR for Landscape Ecology SDM project

#### setup ####
#not a CRAN package, installed from github with devtools::install_github('mrke/NicheMapR')
 # you may need to run R as an admin for that!
 # I also needed to install package 'futile logger'
 # also install global monthly climate means (this is used for the solar radiation calculation if nothing else))
 #get.global.climate(folder="C:/Users/Jordan/Desktop/Landscape ecology/Term_project/Global_Climate")
 # also requires futile.logger package which is on CRAN so can be installed normally



# directories
#sp_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/Species_data/"
#gis_path <- "C:/Users/Jordan/Desktop/Landscape ecology/Term_project/GIS/elevation/"
sp_path <- "E:/LandscapeEco/Species_data/"
gis_path <- "E:/LandscapeEco/GIS/"
out_path <- "E:/LandscapeEco/NicheMapR/"


#packages
library(NicheMapR)
library(lubridate) # deals with dates
library(ggplot2)
library(sp)
library(raster)
library(dismo) # for calculating bioclim variables

# import file with pink lady slipper and background points
pls_spatial <- read.csv(paste(sp_path,"all_pls_points.csv",sep=""))
pls_spatial$X <- NULL
pls_spatial <- SpatialPointsDataFrame(coords=pls_spatial[,c("decimalLongitude","decimalLatitude")],
                                      proj4string=crs("+proj=longlat"),data=pls_spatial)

#import DEM and calculated spatial variables
DEM <- raster(paste(gis_path,"Full_DEM/NASADEM_NC.001_NASADEM_HGT_doy2000042_aid0001.tif",sep=""))
slope <- raster(paste(gis_path,"NYS_slope_degrees",sep=""))
aspect <- raster(paste(gis_path,"NYS_aspect_degrees",sep=""))

# extract values for points
pls_tr <- spTransform(pls_spatial, crs(DEM)) #transform points to crs of DEM

pls_tr$elev <- extract(DEM,pls_tr)
pls_tr$slope <- extract(slope,pls_tr)
pls_tr$aspect <- extract(aspect,pls_tr)

ptdata <- data.frame(pls_tr@data[,c("decimalLongitude","decimalLatitude","elev","slope","aspect","pls_presence")])

#set model parameters
runshade <- 0 #run for only one shade value
runmoist <- 1 #run soil moisture model
soilgrids <- 1 #get soil properties from soilgrids.org
# this increases time from ~18 sec/site to ~32 sec/site
ERR <- 2 # default 1.5 leads to some sites failing



# function to extract data of interest from model result
 
extract_bioclim <- function(micro_out) {
  metdf <- data.frame(micro_out$metout)
  dates <- micro_out$dates
     
  alldf_6h <- cbind(metdf,dates)
  alldf_6h$date <- date(alldf_6h$dates)
  maxt <- aggregate(TALOC ~ date, alldf_6h, max)
  names(maxt) <- c("date","maxt")
  mint <- aggregate(TALOC ~ date, alldf_6h, min)
  names(mint) <- c("date","mint")
  prec <- micro_out$RAINFALL
  
  daily_summary <- merge(maxt,mint)
  daily_summary$prec <- prec
  
  return(daily_summary)
  
}

# things to add 
  # ZH - heat transfer roughness height (m) = 0.02*canopy height
  # D0 - zero plane displacement correction (m) = 0.6* canopy height
  # hori - horizon angles, list of 24
  # look into terrain =1
  # shading - use fPAR?

start_time <- Sys.time()
dem_pt_list <- list()


for(j in 1:length(ptdata$decimalLongitude)) {
  skip <- FALSE
  
  dat <- tryCatch(micro_usa(loc=as.numeric(ptdata[j,c("decimalLongitude","decimalLatitude")]), dstart="01/01/1980", dfinish="31/12/1999",
                            slope=ptdata$slope[j],elev=ptdata$elev[j],aspect=ptdata$aspect[j], #VIEWF=ptdata$svf[j], should be able to incorp but not working
                            runshade=runshade,runmoist=runmoist,soilgrids=soilgrids),
                  error=function(e){skip <<- TRUE})
  
  if(skip==F) {
    dem_pt_list[[j]] <- extract_bioclim(dat) 
    
  } else {
    dem_pt_list[[j]] <- data.frame(date=NA, maxt=NA, mint=NA, prec=NA)
  }
  
  dem_pt_list[[j]]$lat <- as.numeric(ptdata[j,2])
  dem_pt_list[[j]]$lon <- as.numeric(ptdata[j,1])
  dem_pt_list[[j]]$ptID <- j 
  dem_pt_list[[j]]$pls_presence <- ptdata$pls_presence[j]
  
  write.csv(dem_pt_list[[j]],paste(out_path,"microusa_20y_",j,".csv",sep=""))
}

dem_pt_df <- do.call(rbind,dem_pt_list)
write.csv(dem_pt_df,paste(out_path,"20y_Simple_microclimate_summary.csv",sep=""))

end_time <- Sys.time()

start_time
end_time

# summarize by month and calculate bioclim variables
dem_pt_df$month <- month(dem_pt_df$date)
dem_pt_df$year <- year(dem_pt_df$date)

monthly_mean <- aggregate(cbind(maxt,mint)~lat+lon+ptID+pls_presence+month+year,dem_pt_df,mean)
monthly_precip <- aggregate(prec~lat+lon+ptID+pls_presence+month+year,dem_pt_df,sum)

month_dat <- merge(monthly_mean,monthly_precip)
month_dat <- month_dat[order(month_dat$month),] #biovars needs data in order

pt_by_year <- split(month_dat,list(month_dat$ptID,month_dat$year))

bioclim.vars <- data.frame(matrix(ncol=24,nrow=length(pt_by_year)))
names(bioclim.vars) <- c(paste("bio",1:19,sep=""),"lat","lon","ptID","pls_presence","year")

for(i in 1:length(pt_by_year)){
  dat <- pt_by_year[[i]]
  bioclim.vars[i,1:19] <- biovars(prec=dat$prec, tmax=dat$maxt, tmin=dat$mint)
  bioclim.vars[i,]$lat <- dat$lat[1]
  bioclim.vars[i,]$lon <- dat$lon[1]
  bioclim.vars[i,]$ptID <- dat$ptID[1]
  bioclim.vars[i,]$pls_presence <- dat$pls_presence[1]
  bioclim.vars[i,]$year <- dat$year[1]
}

bioclim.vars$year <- NULL
mean_bioclim <- aggregate(.~lat+lon+ptID+pls_presence,bioclim.vars,mean)

write.csv(mean_bioclim,paste(sp_path,"bioclim_20ymean_simplemicroclimate.csv",sep=""))
