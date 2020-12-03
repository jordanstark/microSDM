# figure comparing bioclim and nichemapr predictions


library(ggplot2)
library(lubridate) #dates
library(sp) #GIS (spatial points)
library(raster) #GIS (rasters)
library(tmap) # for plotting maps
library(tmaptools) #also maps
library(tidyr) #mostly for the pivot_longer tranformation function
library(stringr) # dealing with text
library(patchwork) # for plotting multiple figures
library(car) # for anova
library(rstatix) # wrapper for car for repeated measures

states <- getData("GADM",country="USA",level=1) #download state data
NY <- states[states$NAME_1=="New York",]
NY_poly <- SpatialPolygons(NY@polygons,proj4string=crs("+proj=longlat"))


micro <- read.csv("c:/Users/Jordan/Desktop/Landscape ecology/Term_project/bioclim_20ymean_shademicroclimate.csv")
micro$X <- NULL
micro$pt_type <- NA
micro$pt_type[micro$pls_presence == F] <- "bg"
micro$pt_type[micro$pls_presence == T] <- "pls"
micro$pls_presence <- NULL
micro$ptID <- paste(micro$ptID,"pls",sep="")

micro_rbs <- read.csv("c:/Users/Jordan/Desktop/Landscape ecology/Term_project/bioclim_20ymean_shademicroclimate_rbs.csv")
micro_rbs$X <- NULL
micro_rbs$pt_type <- NA
micro_rbs$pt_type[micro_rbs$rbs_presence == F] <- "bg"
micro_rbs$pt_type[micro_rbs$rbs_presence == T] <- "rbs"
micro_rbs$rbs_presence <- NULL
micro_rbs$ptID <- paste(micro_rbs$ptID,"rbs",sep="")

all_micro <- rbind(micro,micro_rbs)

micro_sp <- SpatialPointsDataFrame(coords=cbind(all_micro$lon,all_micro$lat),
                                   proj4string=crs("+proj=longlat"),data=all_micro)

bioclim_30s <- getData('worldclim',var="bio",res=0.5,lat=42,lon=-75,)
bioclim_2.5m <- getData('worldclim',var="bio",res=2.5,lat=42,lon=-75,)


bioclim_30s_sp <- data.frame(raster::extract(bioclim_30s,micro_sp))
names(bioclim_30s_sp) <- paste("bio",1:19,sep="")
bioclim_30s_sp$lat <- micro_sp$lat
bioclim_30s_sp$lon <- micro_sp$lon
bioclim_30s_sp$ptID <- micro_sp$ptID
bioclim_30s_sp$pt_type <- micro_sp$pt_type

bioclim_2.5m_sp <- data.frame(raster::extract(bioclim_2.5m,micro_sp))
names(bioclim_2.5m_sp) <- paste("bio",1:19,sep="")
bioclim_2.5m_sp$lat <- micro_sp$lat
bioclim_2.5m_sp$lon <- micro_sp$lon
bioclim_2.5m_sp$ptID <- micro_sp$ptID
bioclim_2.5m_sp$pt_type <- micro_sp$pt_type


# temperature values are *10, correcting
bioclim_30s_sp$bio1 <- bioclim_30s_sp$bio1/10
bioclim_30s_sp$bio2 <- bioclim_30s_sp$bio2/10
bioclim_30s_sp$bio3 <- bioclim_30s_sp$bio3/100 #bio 3 and 4 are x100
bioclim_30s_sp$bio4 <- bioclim_30s_sp$bio4/100
bioclim_30s_sp$bio5 <- bioclim_30s_sp$bio5/10
bioclim_30s_sp$bio6 <- bioclim_30s_sp$bio6/10
bioclim_30s_sp$bio7 <- bioclim_30s_sp$bio7/10
bioclim_30s_sp$bio8 <- bioclim_30s_sp$bio8/10
bioclim_30s_sp$bio9 <- bioclim_30s_sp$bio9/10
bioclim_30s_sp$bio10 <- bioclim_30s_sp$bio10/10

bioclim_2.5m_sp$bio1 <- bioclim_2.5m_sp$bio1/10
bioclim_2.5m_sp$bio2 <- bioclim_2.5m_sp$bio2/10
bioclim_2.5m_sp$bio3 <- bioclim_2.5m_sp$bio3/100 #bio 3 and 4 are x100
bioclim_2.5m_sp$bio4 <- bioclim_2.5m_sp$bio4/100
bioclim_2.5m_sp$bio5 <- bioclim_2.5m_sp$bio5/10
bioclim_2.5m_sp$bio6 <- bioclim_2.5m_sp$bio6/10
bioclim_2.5m_sp$bio7 <- bioclim_2.5m_sp$bio7/10
bioclim_2.5m_sp$bio8 <- bioclim_2.5m_sp$bio8/10
bioclim_2.5m_sp$bio9 <- bioclim_2.5m_sp$bio9/10
bioclim_2.5m_sp$bio10 <- bioclim_2.5m_sp$bio10/10



micro_long <- pivot_longer(micro_sp@data,cols=starts_with("bio"),values_to="micro_value",names_to="biovar")
bioclim_30s_long <- pivot_longer(bioclim_30s_sp,cols=starts_with("bio"),values_to="bio_30s_value",names_to="biovar")
bioclim_2.5m_long <- pivot_longer(bioclim_2.5m_sp,cols=starts_with("bio"),values_to="bio_2.5m_value",names_to="biovar")


all_dat <- merge(micro_long,bioclim_30s_long)
all_dat <- merge(all_dat,bioclim_2.5m_long)

#all_dat$micro_delta <- abs(all_dat$micro_value - all_dat$bio_value)

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

# ggplot(all_dat, aes(x=pt_type, y=micro_delta)) +
#   geom_boxplot() +
#   facet_wrap(~biovar,scales="free") +
#   theme_classic()

all_dat$frac.delta.30s <- (all_dat$bio_30s_value - all_dat$bio_2.5m_value) / all_dat$bio_2.5m_value
all_dat$frac.delta.micro <- (all_dat$micro_value - all_dat$bio_2.5m_value) / all_dat$bio_2.5m_value

ggplot(all_dat,aes(x=pt_type,y=frac.delta.micro,color=pt_type)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  facet_wrap(~biovar,scales="free",nrow=2) +
  theme_classic()

longer_dat <- pivot_longer(all_dat,cols=c("micro_value","bio_30s_value","bio_2.5m_value"),
                           names_to="scale",values_to="value")

longer_dat <- data.frame(longer_dat)
longer_dat$scale[longer_dat$scale=="micro_value"] <- "micro"
longer_dat$scale[longer_dat$scale=="bio_30s_value"] <- "macro"
longer_dat$scale[longer_dat$scale=="bio_2.5m_value"] <- "regional"

longer_dat$scale <- factor(longer_dat$scale,levels=c("regional","macro","micro"),ordered=T)


se <- function(x) {mean(x,na.rm=T)/sum(!is.na(x))}

byvar <- split(longer_dat,longer_dat$biovar)
plots <- list()
stat.results <- data.frame(matrix(nrow=19,ncol=7))
names(stat.results) <- c("biovar","pls.pt_type","pls.scale","pls.intxn","rbs.pt_type","rbs.scale","rbs.intxn")

for(i in 1:length(byvar)){
  dat <- byvar[[i]]
  dat_pls <- dat[dat$pt_type != "rbs",]
  dat_rbs <- dat[dat$pt_type != "pls",]
  
  pls.test <- anova_test(value ~ pt_type*scale,
                         data=dat_pls,
                         dv=value, #dep var = climate values
                         wid=ptID, #repeated measures ID
                         type=3)
  
  rbs.test <- anova_test(value ~ pt_type*scale,
                          data=dat_rbs,
                          dv=value, #dep var = climate values
                          wid=ptID, #repeated measures ID
                          type=3) 
  
  stat.results[i,1] <- dat$biovar[1]
  stat.results[i,2:4] <- pls.test$p
  stat.results[i,5:7] <- rbs.test$p
  
  dat_mean <- aggregate(value ~ scale + pt_type,dat,mean)
  names(dat_mean) <- c("scale","pt_type","mean_value")
  dat_se <- aggregate(value ~ scale + pt_type,dat,se)
  names(dat_se) <- c("scale","pt_type","se_value")
  
  dat_summary <- merge(dat_mean,dat_se)
  
  plots[[i]]  <-     ggplot(dat_summary,aes(x=scale,color=pt_type)) +
                        #geom_line(data=dat,aes(x=scale,y=value,color=pt_type,group=ptID),alpha=0.05,inherit.aes=F) +
                        #geom_violin(data=dat,aes(x=scale,y=value,color=pt_type,fill=pt_type),position="identity",alpha=0.1) +
                        geom_errorbar(aes(ymin=mean_value-se_value,ymax=mean_value+se_value),width=0.2,size=1) +
                        geom_line(aes(y=mean_value,group=pt_type),size=1) +
                        scale_color_manual(values=c("#666699","#33CC33","#CC0000")) +
                        theme_classic() +
                        xlab("") + ylab(paste(dat$text_name[1])) +
                        theme(legend.position="none")
                        


  }


stat.results

stat.sig <- data.frame(stat.results[,2:7] < (0.05/19)) #bonferroni correction for 19 tests
stat.sig$biovar <- stat.results[,1]
stat.sig <- merge(stat.sig,bioclim_names)

sig.pls.intxn <- stat.sig$biovar[stat.sig$pls.intxn == T]
sig.rbs.intxn <- stat.sig$biovar[stat.sig$rbs.intxn == T]

sig.plots <- plots[ which(stat.sig$pls.intxn==T | stat.sig$rbs.intxn == T)]

plots <- list()
iter <- 1

for(i in which(stat.sig$pls.intxn == T)) {
  pnum <- iter
  dat <- byvar[[i]]
  dat <- dat[dat$pt_type != "rbs",]
  
  dat_mean <- aggregate(value ~ scale + pt_type,dat,mean)
  names(dat_mean) <- c("scale","pt_type","mean_value")
  dat_se <- aggregate(value ~ scale + pt_type,dat,se)
  names(dat_se) <- c("scale","pt_type","se_value")
  
  dat_summary <- merge(dat_mean,dat_se)
  
  plots[[pnum]]  <-   ggplot(dat_summary,aes(x=scale,color=pt_type)) +
                          geom_line(data=dat,aes(x=scale,y=value,color=pt_type,group=ptID),alpha=0.05,inherit.aes=F) +
                          #geom_violin(data=dat,aes(x=scale,y=value,color=pt_type,fill=pt_type),position="identity",alpha=0.1) +
                          geom_errorbar(aes(ymin=mean_value-se_value,ymax=mean_value+se_value),width=0.2,size=1) +
                          geom_line(aes(y=mean_value,group=pt_type),size=1) +
                          scale_color_manual(values=c("black","deeppink")) +
                          theme_classic() +
                          xlab("") + ylab(paste(dat$text_name[1])) +
                          theme(legend.position="none")
  
  iter <- iter+1
                            

}

rbs.plots <- list()
for(i in which(stat.sig$rbs.intxn == T)) {
  pnum <- iter
  dat <- byvar[[i]]
  dat <- dat[dat$pt_type != "pls",]
  
  dat_mean <- aggregate(value ~ scale + pt_type,dat,mean)
  names(dat_mean) <- c("scale","pt_type","mean_value")
  dat_se <- aggregate(value ~ scale + pt_type,dat,se)
  names(dat_se) <- c("scale","pt_type","se_value")
  
  dat_summary <- merge(dat_mean,dat_se)
  
  plots[[pnum]]  <-   ggplot(dat_summary,aes(x=scale,color=pt_type)) +
    geom_line(data=dat,aes(x=scale,y=value,color=pt_type,group=ptID),alpha=0.05,inherit.aes=F) +
    #geom_violin(data=dat,aes(x=scale,y=value,color=pt_type,fill=pt_type),position="identity",alpha=0.1) +
    geom_errorbar(aes(ymin=mean_value-se_value,ymax=mean_value+se_value),width=0.2,size=1) +
    geom_line(aes(y=mean_value,group=pt_type),size=1) +
    scale_color_manual(values=c("black","red2")) +
    theme_classic() +
    xlab("") + ylab(paste(dat$text_name[1])) +
    theme(legend.position="none")
  
  iter <- iter+1
  
  
}

all_plots <- list(pls.plots,rbs.plots)

wrap_plots(plots,nrow=2)

summary_dat <- aggregate(cbind(micro_value,bio_30s_value,bio_2.5m_value,frac.delta.30s,frac.delta.micro) ~ biovar + pt_type, 
                         all_dat, mean)



summary_long <- pivot_longer(summary_dat,cols=c("micro_value","bio_30s_value","bio_2.5m_value"),
                             values_to="bio_val",names_to="climate_mod")
summary_long$climate_mod <- factor(summary_long$climate_mod)


ggplot(summary_long,aes(x=climate_mod,y=bio_val,color=pt_type)) +
  geom_point() +
  geom_line() +
  facet_wrap(~biovar,scales="free") +
  theme_classic()


var_split <- split(all_dat,all_dat$biovar)



bio1 <- var_split[[1]]

bio1$frac.delta.30s <- (bio1$bio_30s_value - bio1$bio_2.5m_value) / bio1$bio_2.5m_value
bio1$frac.delta.micro <- (bio1$micro_value - bio1$bio_2.5m_value) / bio1$bio_2.5m_value

ggplot(bio1,aes(x=pls_presence, y=frac.delta.micro, color=pls_presence)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) +
  theme_classic()

bio1_long <- pivot_longer(bio1,cols=c("micro_value","bio_30s_value","bio_2.5m_value"),
                          values_to="bio1",names_to="climate_mod")

ggplot(bio1_long,aes(x=climate_mod,y=bio1)) +
  geom_boxplot(aes(color=pls_presence),outlier.shape=NA) +
  geom_jitter(aes(color=pls_presence),width=0.1) +
  #geom_line(aes(group=ptID,color=pls_presence), alpha=0.5) +
  theme_classic()



############################################ not updated #########################

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
