##Author: Mike Back
##Organization: Kent State University
##Project: FoSTER Project, 2023 manuscript
##Script purpose: Near surface bulk density and Ksat correlation stats
##Date updated:7/26/2023


##Import data
#call "near_surface_bulk_density.csv" for bulk density data
#call "Ksat.csv" for Ksat data
density <- read.csv("C:\\Users\\mikeb\\Documents\\github\\foster_soil\\data\\near_surface_bulk_density.csv", stringsAsFactors = T)
Ksat <- read.csv("C:\\Users\\mikeb\\Documents\\github\\foster_soil\\data\\ksat.csv", stringsAsFactors = T)

ksat.prerip <- subset(Ksat, pre_post=="pre"|pre_post=="no")#no = NA (i.e. reference sites or mine sites that haven't been ripped yet)
density.prerip <- subset(density, rip=="pre"|rip=="no")

merged <- merge(density.prerip, ksat.prerip, by="sample")
View(merged)
summary(merged)

#transform ksat
merged$log.ksat <- log(merged$Ksat.ms)

#check normality
shapiro.test(merged$log.ksat)$p
shapiro.test(merged$bulk_density)$p
hist(merged$log.ksat)#how do I account for the 3 truncated measurements of 1.47x10-8?

#regression model
density.ksat.lm <- lm(log.ksat~bulk_density, data = merged)
summary(density.ksat.lm)
plot(density.ksat.lm, which=c(1,2))
plot(log.ksat~bulk_density, data = merged)

#correlation
cor.test(merged$Ksat.ms, merged$bulk_density)#without log transformed Ksat
cor.test(merged$log.ksat, merged$bulk_density)

#at DV only 
merged.DV <- subset(merged,  site.x=="DV")
merged.DV$log.ksat <- log(merged.DV$Ksat.ms)
View(merged.DV)

DV.lm <- lm(log.ksat~bulk_density, data = merged.DV)
summary(DV.lm)

cor.test(merged.DV$log.ksat, merged.DV$bulk_density)

#at SQ only
merged.SQ <- subset(merged,  site.x=="SQ")
merged.SQ$log.ksat <- log(merged.SQ$Ksat.ms)
View(merged.SQ)

SQ.lm <- lm(log.ksat~bulk_density, data = merged.SQ)
summary(SQ.lm)

cor.test(merged.SQ$log.ksat, merged.SQ$bulk_density)
