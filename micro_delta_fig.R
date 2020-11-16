# figure comparing bioclim and nichemapr predictions


library(ggplot2)
library(lubridate)
library(sp)
library(raster)
library(tmap)
library(tmaptools)
library(tidyr)
library(stringr)
library(patchwork)

states <- getData("GADM",country="USA",level=1) #download state data
NY <- states[states$NAME_1=="New York",]
NY_poly <- SpatialPolygons(NY@polygons,proj4string=crs("+proj=longlat"))


micro <- read.csv("c:/Users/Jordan/Desktop/Landscape ecology/Term_project/bioclim_20ymean_shademicroclimate.csv")
micro$X <- NA

micro_sp <- SpatialPointsDataFrame(coords=cbind(micro$lon,micro$lat),proj4string=crs("+proj=longlat"),data=micro)

bioclim_1km <- getData('worldclim',var="bio",res=0.5,lat=42,lon=-75,)
#mint_1km <- getData('worldclim',var="tmin",res=0.5,lat=-42,lon=-75)

bioclim_sp <- data.frame(raster::extract(bioclim_1km,micro_sp))
names(bioclim_sp) <- paste("bio",1:19,sep="")
bioclim_sp$lat <- micro_sp$lat
bioclim_sp$lon <- micro_sp$lon
bioclim_sp$ptID <- micro_sp$ptID
bioclim_sp$pls_presence <- micro_sp$pls_presence

# temperature values are *10, correcting
bioclim_sp$bio1 <- bioclim_sp$bio1/10
bioclim_sp$bio2 <- bioclim_sp$bio2/10
bioclim_sp$bio3 <- bioclim_sp$bio3/100 #bio 3 and 4 are x100
bioclim_sp$bio4 <- bioclim_sp$bio4/100
bioclim_sp$bio5 <- bioclim_sp$bio5/10
bioclim_sp$bio6 <- bioclim_sp$bio6/10
bioclim_sp$bio7 <- bioclim_sp$bio7/10
bioclim_sp$bio8 <- bioclim_sp$bio8/10
bioclim_sp$bio9 <- bioclim_sp$bio9/10
bioclim_sp$bio10 <- bioclim_sp$bio10/10



micro_long <- pivot_longer(micro_sp@data,cols=starts_with("bio"),values_to="micro_value",names_to="biovar")
bioclim_long <- pivot_longer(bioclim_sp,cols=starts_with("bio"),values_to="bio_value",names_to="biovar")

all_dat <- merge(micro_long,bioclim_long)
all_dat$micro_delta <- abs(all_dat$micro_value - all_dat$bio_value)

bioclim_names <- c(bio1="Annual Temp",
                   bio2="Diurnal Temp Range",
                   bio3="Isothermality",
                   bio4="Temp Seasonality",
                   bio5="Max Temp",
                   bio6="Min Temp",
                   bio7="Annual Temp Range",
                   bio8="Mean Temp Wettest Quarter",
                   bio9="Mean Temp Driest Quarter",
                   bio10="Mean Temp Warmest Quarter",
                   bio11="Mean Temp Coldest Quarter",
                   bio12="Annual Precip",
                   bio13="Precip Wettest Month",
                   bio14="Precip Driest Month",
                   bio15="Precip Seasonality",
                   bio16="Precip Wettest Quarter",
                   bio17="Precip Driest Quarter",
                   bio18="Precip Warmest Quarter",
                   bio19="Precip Coldest Quarter")
bioclim_names <- data.frame(biovar=names(bioclim_names),text_name=bioclim_names)

all_dat <- merge(all_dat,bioclim_names)

ggplot(all_dat, aes(x=pls_presence, y=micro_delta)) +
  geom_boxplot() +
  facet_wrap(~biovar,scales="free") +
  theme_classic()


var_split <- split(all_dat,all_dat$biovar)

tests <- data.frame(biovar=NA, text_name=NA,
                    mean_bioclim_pls=NA,mean_bioclim_bg=NA,
                    mean_microclim_pls=NA,mean_microclim_bg=NA,
                    sd_bioclim_pls=NA,sd_bioclim_bg=NA,
                    sd_microclim_pls=NA,sd_microclim_bg=NA,
                    p_diff=NA,p_pls=NA,
                    p_bioclim=NA,p_microclim=NA)

plot_plseff <- list()
plot_climdiff <- list()
plot_bios <- list()
plot_micros <- list()


for(i in 1:length(var_split)) {
  dat <- var_split[[i]]
  biovar <- dat$biovar[1]
  text_name <- dat$text_name[1]
  mean_bioclim_pls <- mean(dat$bio_value[dat$pls_presence==T],na.rm=T)
  mean_bioclim_bg <- mean(dat$bio_value[dat$pls_presence==F],na.rm=T)
  mean_microclim_pls <- mean(dat$micro_value[dat$pls_presence==T],na.rm=T)
  mean_microclim_bg <- mean(dat$micro_value[dat$pls_presence==F],na.rm=T)
  
  sd_bioclim_pls <- sd(dat$bio_value[dat$pls_presence==T],na.rm=T)
  sd_bioclim_bg <- sd(dat$bio_value[dat$pls_presence==F],na.rm=T)
  sd_microclim_pls <- sd(dat$micro_value[dat$pls_presence==T],na.rm=T)
  sd_microclim_bg <- sd(dat$micro_value[dat$pls_presence==F],na.rm=T)
  
  diff.test <- t.test(x=dat$micro_value,y=dat$bio_value,paired=T)
  p_diff <- round(diff.test$p.value,digits=5)
  
  sp.test <- t.test(dat$micro_delta ~ dat$pls_presence)
  p_pls <- round(sp.test$p.value,digits=5)
  
  bioclim.test <- t.test(dat$bio_value ~ dat$pls_presence)
  p_bioclim <- round(bioclim.test$p.value,digits=5)
  
  microclim.test <- t.test(dat$micro_value ~ dat$pls_presence)
  p_microclim <- round(microclim.test$p.value,digits=5)
  
  tests[i,"biovar"] <- biovar
  tests[i,"text_name"] <- text_name
  
  tests[i,3:14] <- c(mean_bioclim_pls,mean_bioclim_bg,
                 mean_microclim_pls,mean_microclim_bg,
                 sd_bioclim_pls,sd_bioclim_bg,
                 sd_microclim_pls,sd_microclim_bg,
                 p_diff,p_pls,p_bioclim,p_microclim)
  
  plot_micros[[i]] <- ggplot(dat,aes(y=micro_value,x=pls_presence)) +
                        geom_boxplot(outlier.shape=NA) +
                        geom_jitter(aes(color=pls_presence),width=0.2) +
                        ggtitle(paste(text_name,"     ","p=",p_microclim,sep="")) +
                        theme_classic() +
                        theme(legend.position="none") +
                        xlab("") + ylab("")
  
  plot_bios[[i]] <- ggplot(dat,aes(y=bio_value,x=pls_presence)) +
                            geom_boxplot(outlier.shape=NA) +
                            geom_jitter(aes(color=pls_presence),width=0.2) +
                            ggtitle(paste(text_name,"     ","p=",p_bioclim,sep="")) +
                            theme_classic() +
                            theme(legend.position="none") +
                            xlab("") + ylab("")
                      
  
  plot_plseff[[i]] <-   ggplot(dat,aes(y=micro_delta,x=pls_presence,group=pls_presence)) +
                          geom_boxplot(outlier.shape=NA) +
                          geom_jitter(aes(color=pls_presence),width=0.2) +
                          ggtitle(paste(text_name,"     ","p=",p_pls,sep="")) +
                          theme_classic() +
                          theme(legend.position="none") +
                          xlab("") + ylab("")
  
  dat_long <- pivot_longer(dat,cols=c("micro_value","bio_value"),names_to="clim",values_to="value")
  
  plot_climdiff[[i]] <- ggplot(dat_long,aes(x=clim,y=value)) +
                          geom_point(aes(color=pls_presence)) +
                          geom_line(aes(group=ptID,color=pls_presence),alpha=0.5) +
                          ggtitle(paste(text_name)) +
                          theme_classic() +
                          theme(legend.position="none") +
                          xlab("") + ylab("")
  

}


sig_plots <- plot_plseff[which(tests$p_pls<(0.05/19))] # bonferroni correction

wrap_plots(sig_plots,nrow=1)

sig_bios <- plot_bios[which(tests$p_bioclim<(0.05/19))]
wrap_plots(sig_bios,nrow=1)

sig_micros <- plot_micros[which(tests$p_microclim<(0.05/19))]
wrap_plots(sig_micros,nrow=1)



impt_vars <- c("bio1","bio5","bio6","bio12","bio13","bio14")

good_plots <- plot_climdiff[which(tests$biovar %in% impt_vars)]

wrap_plots(good_plots,nrow=1)


sig_pls <- tests[tests$p_pls < 0.05,] 

sig_pls_long <- pivot_longer(sig_pls,cols=2:9,names_to="var",values_to="val")

sig_pls_long$clim <- str_split(sig_pls_long$var,pattern="_",simplify=T)[,2]  
sig_pls_long$pls <- str_split(sig_pls_long$var,pattern="_",simplify=T)[,3]  
sig_pls_long$type <- str_split(sig_pls_long$var,pattern="_",simplify=T)[,1]  




ggplot(sig_pls,aes(x=biovar)) +
  geom_point(aes(y=mean_bioclim_pls),color="blue") +
  geom_point(aes(y=mean_bioclim_bg),color="red") +
  geom_point(aes(y=mean_microclim_pls),color="lightblue") +
  geom_point(aes(y=mean_microclim_bg),color="pink")
