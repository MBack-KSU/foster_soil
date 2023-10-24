##Author: Mike Back
##Organization: Kent State University
##Project: FoSTER Project, 2023 manuscript
##Script purpose: Spatial stats and semivariogram figures
##Date updated:7/26/2023


cvnp.data <- read.csv("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\spatial_and_bulk_density_data.csv")

#subset mine sites for pre-rip bulk density omitting NA's and repeats (which occur in post-rip and cause issues)
DV.mine <- subset(cvnp.data, site=="DV")

SQ <- subset(cvnp.data, site=="SQ")
SQ.no.na <- na.omit(SQ[,1:6])
SQ.mine <- SQ.no.na[c(1:10,12:21),]

HH.mine <- subset(cvnp.data, site=="HH")
HH.mine <- na.omit(HH.mine[,1:6])

CT.mine <- subset(cvnp.data, site=="CT")
CT.mine <- na.omit(CT.mine[,1:6])

RR.mine <- subset(cvnp.data, site=="RR")
RR.mine <- na.omit(RR.mine[,1:6])

#####
##Run Moran's I for spatial autocorrelation at each mine site for pre-rip bulk denisty
library(ape)#use this for global moran's I

#DV Moran.I test using ape package (global)
DV.dists <- as.matrix(dist(cbind(DV.mine$longitude, DV.mine$latitude))) #measure distances based on coords
DV.dists.inv <- 1/DV.dists #take inverse of the distances
diag(DV.dists.inv) <- 0 #set diagonal values in matrix to 0
DV.dists.inv[1:5, 1:5] 

Moran.I(DV.mine$pre_rip_bulk_density, DV.dists.inv)#function in ape for global morans I

#Moran's I test using nearest neighbor
library(spdep)#use this for nearest neighbors
#this chunk of code converts the DV.mine df to a SpatialPointsDataFrame with utm coordinates
coordinates(DV.mine) <- ~longitude+latitude
proj4string(DV.mine) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
DV.utm <- spTransform(DV.mine, CRS("+proj=utm +zone=17 ellps=WGS84"))
class(DV.mine)
DV.utm

DV.nb <- dnearneigh(coordinates(DV.utm), 0, 350)#0 is min. dist, and the second number is max. dist between neighbors
DV.weights <- nb2listwdist(DV.nb, DV.mine)
moran.test(DV.mine$pre_rip_bulk_density,DV.weights)

#SQ Moran.I test
SQ.dists <- as.matrix(dist(cbind(SQ.mine$longitude, SQ.mine$latitude))) #measure distances based on coords
SQ.dists.inv <- 1/SQ.dists #take invert of the distances
diag(SQ.dists.inv) <- 0 #set the inverse to 0
SQ.dists.inv[1:5, 1:5]

Moran.I(SQ.mine$pre_rip_bulk_density, SQ.dists.inv)

coordinates(SQ.mine) <- ~longitude+latitude # this comes from the {sp} package
proj4string(SQ.mine) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
SQ.utm <- spTransform(SQ.mine, CRS("+proj=utm +zone=17 ellps=WGS84"))
class(SQ.mine)
SQ.utm

SQ.nb <- dnearneigh(coordinates(SQ.utm), 0, 350)
SQ.weights <- nb2listwdist(SQ.nb, SQ.mine)
moran.test(SQ.mine$pre_rip_bulk_density,SQ.weights)

#HH Moran's test
HH.dists <- as.matrix(dist(cbind(HH.mine$longitude, HH.mine$latitude))) #measure distances based on coords
HH.dists.inv <- 1/HH.dists #take invert of the distances
diag(HH.dists.inv) <- 0 #set the inverse to 0
HH.dists.inv[1:5, 1:5]

Moran.I(HH.mine$pre_rip_bulk_density, HH.dists.inv)

coordinates(HH.mine) <- ~longitude+latitude # this comes from the {sp} package
proj4string(HH.mine) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
HH.utm <- spTransform(HH.mine, CRS("+proj=utm +zone=17 ellps=WGS84"))
class(HH.mine)
HH.utm

HH.nb <- dnearneigh(coordinates(HH.utm), 0, 350)
HH.weights <- nb2listwdist(HH.nb, HH.mine)
moran.test(HH.mine$pre_rip_bulk_density,HH.weights)

#CT Moran's test
CT.dists <- as.matrix(dist(cbind(CT.mine$longitude, CT.mine$latitude))) #measure distances based on coords
CT.dists.inv <- 1/CT.dists #take invert of the distances
diag(CT.dists.inv) <- 0 #set the inverse to 0
CT.dists.inv[1:5, 1:5]

Moran.I(CT.mine$pre_rip_bulk_density, CT.dists.inv)

coordinates(CT.mine) <- ~longitude+latitude # this comes from the {sp} package
proj4string(CT.mine) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
CT.utm <- spTransform(CT.mine, CRS("+proj=utm +zone=17 ellps=WGS84"))
class(CT.mine)
CT.utm

CT.nb <- dnearneigh(coordinates(CT.utm), 0, 350)
CT.weights <- nb2listwdist(CT.nb, CT.mine)
moran.test(CT.mine$pre_rip_bulk_density,CT.weights)

#RR moran's test
RR.dists <- as.matrix(dist(cbind(RR.mine$longitude, RR.mine$latitude))) #measure distances based on coords
RR.dists.inv <- 1/RR.dists #take invert of the distances
diag(RR.dists.inv) <- 0 #set the inverse to 0
RR.dists.inv[1:5, 1:5]

Moran.I(RR.mine$pre_rip_bulk_density, RR.dists.inv)

coordinates(RR.mine) <- ~longitude+latitude # this comes from the {sp} package
proj4string(RR.mine) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
RR.utm <- spTransform(RR.mine, CRS("+proj=utm +zone=17 ellps=WGS84"))
class(RR.mine)
RR.utm

RR.nb <- dnearneigh(coordinates(RR.utm), 0, 350)
RR.weights <- nb2listwdist(RR.nb, RR.mine)
moran.test(RR.mine$pre_rip_bulk_density,RR.weights)

#####
#Make semivariograms of pre-rip bulk density for each site
library(rgdal)
library(gstat)
library(sp)
library(vegan)
library(MASS)
library(ggplot2)

##DV variogram
DV.vec=as.vector(c(40,55,80,85,107,115,130,150,160,189,215,250,350))#Set up bins, tried to equally distribute pairs across bins
DV.vgm=variogram(pre_rip_bulk_density~1,DV.utm,boundaries=DV.vec)
DV.vgm

DV.plot <- ggplot(DV.vgm, aes(x = dist, y = gamma)) + 
  geom_point(size=2) + 
  xlim(0,300) + 
  ylim(0,0.06) + 
  xlab("Distance (m)") + 
  ylab("Semivariance") + 
  ggtitle("DV") + 
  theme_classic() + 
  theme(text = element_text(size = 10))
DV.plot
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\semivariograms\\DV_vgm.jpeg", plot=DV.plot, height=5, width=7, units=c("cm"), dpi=600)
# if we want to add moran's I "annotate(geom = "text", x=25, y=0.035, label="Moran's I = -0.058; p = 0.262")"

##SQ variogram
#majority same code as for DV
SQ.vec=as.vector(c(52,65,77,84,105,115,130,150,160,189,215,250,350))
SQ.vgm=variogram(pre_rip_bulk_density~1,SQ.utm,boundaries = SQ.vec)
SQ.vgm

SQ.plot <- ggplot(SQ.vgm, aes(x = dist, y = gamma)) + 
  geom_point(size=2) +
  xlim(0,300) + ylim(0,0.06) + 
  xlab("Distance (m)") + 
  ylab("Semivariance") + 
  ggtitle("SQ") + 
  theme_classic() + 
  theme(text = element_text(size = 10))
SQ.plot
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\semivariograms\\SQ_vgm.jpeg", plot=SQ.plot, height=5, width=7, units=c("cm"), dpi=600)


##CT variogram
#majority same code as for DV
CT.vec=as.vector(c(37,50,67,81,100,117,140,350))
CT.vgm=variogram(pre_rip_bulk_density~1,CT.utm,boundaries = CT.vec)
CT.vgm
summary(CT.vgm)
str(CT.vgm)

CT.plot <- ggplot(CT.vgm, aes(x = dist, y = gamma)) +
  geom_point(size=2) +
  ylim(0,0.15) +
  xlab("Distance (m)") + 
  ylab("Semivariance") + 
  ggtitle("CT") + 
  theme_classic() + 
  theme(text = element_text(size = 10))+
  scale_x_continuous(breaks=seq(0,200, by=100),
                     limits = c(0,200))
CT.plot
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\semivariograms\\CT_vgm.jpeg", plot=CT.plot, height=5, width=5, units=c("cm"), dpi=600)


##HH variogram
#majority same code as for DV
HH.vec=as.vector(c(35,52,72,86,105,130,160,190,215,300))
HH.vgm=variogram(pre_rip_bulk_density~1,HH.utm,boundaries = HH.vec)
HH.vgm

HH.plot <- ggplot(HH.vgm, aes(x = dist, y = gamma)) + 
  geom_point(size=2) +
  xlim(0,300) + ylim(0,0.06) + 
  xlab("Distance (m)") + 
  ylab("Semivariance") + 
  ggtitle("HH") + 
  theme_classic() + 
  theme(text = element_text(size = 10))
HH.plot
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\semivariograms\\HH_vgm.jpeg", plot=HH.plot, height=5, width=7, units=c("cm"), dpi=600)


##RR variogram
#majority same code as for DV
RR.vec=as.vector(c(48,65,85,102,125,150,300))
RR.vgm=variogram(pre_rip_bulk_density~1,RR.utm,boundaries = RR.vec)
RR.vgm

RR.plot <- ggplot(RR.vgm, aes(x = dist, y = gamma)) + 
  geom_point(size=2) +
  ylim(0,0.06) + 
  xlab("Distance (m)") + 
  ylab("Semivariance") + 
  ggtitle("RR") + 
  theme_classic() + 
  theme(text = element_text(size = 10)) + 
  scale_x_continuous(breaks=seq(0,200, by=100),
                     limits = c(0,200))
RR.plot
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\semivariograms\\RR_vgm.jpeg", plot=RR.plot, height=5, width=5, units=c("cm"), dpi=600)
