##Author: Mike Back
##Organization: Kent State University
##Project: FoSTER Project, 2023 manuscript
##Script purpose: Near sufrace bulk density stats and figures
##Date updated:7/26/2023


density <- read.csv("C:\\Users\\mikeb\\Documents\\github\\foster_soil\\data\\near_surface_bulk_density.csv", stringsAsFactors = T)
View(density)
density$number <- factor(density$number)
summary(density)

library("nlme")
library("lme4")
library("lmerTest")
library("emmeans")
library("ggplot2")
library("MuMIn")

##create subsets
CT.RR.CTR <- subset(density, site=="CT"|site=="CTR"|site=="RR")
CT.CTR <- subset(density, site=="CT"|site=="CTR")
CT <- subset(density, site=="CT")
CTR <- subset(density, site=="CTR")
RR.CTR <- subset(density, site=="RR"|site=="CRR")
RR <- subset(density, site=="RR")
HH.HHR <- subset(density, site=="HH"|site=="HHR")
HH <- subset(density, site=="HH")
HHR <- subset(density, site=="HHR")
DVR <- subset(density, site=="DVR")
DV.DVR <- subset(density, site=="DV"|site=="DVR")
SQR <- subset(density, site=="SQR")
SQ.SQR <- subset(density, site=="SQ"|site=="SQR")
#subset for pre vs ref t-test
DVpre.DVR <- subset(DV.DVR, rip=="pre"|site=="DVR")
SQpre.SQR <- subset(SQ.SQR, rip=="pre"|site=="SQR")
##subset for pre-post analysis
DV.SQ.PrePost <- subset(density, site=="DV"|site=="SQ")
summary(SQ.SQR)
##subset for mine vs reference
mine.ref <- subset(density, rip=="no"|rip=="pre")
View(mine.ref)


##test for all mine vs reference
mine.ref.LMER <- lmer(bulk_density~type + (1 | general), data = mine.ref)
plot(mine.ref.LMER)
qqnorm(resid(mine.ref.LMER))
anova(mine.ref.LMER)
summary(mine.ref.LMER)
emm1=emmeans(mine.ref.LMER, ~type)
pairs(emm1, adjust = "holm")
r.squaredGLMM(mine.ref.LMER)
##model without general as random effect
mine.ref.LM <- lm(bulk_density~type, data = mine.ref)
plot(mine.ref.LM)
qqnorm(resid(mine.ref.LM))
anova(mine.ref.LM)
summary(mine.ref.LM)
##likelihood ratio test between the two models
anova(mine.ref.LMER, mine.ref.LM)

##need to find which sites are significantly different for boxplot
mine.ref.sites <- lm(bulk_density~site, data=mine.ref)
anova(mine.ref.sites)
emm3=emmeans(mine.ref.sites, ~site)
pairs(emm3, adjust = "holm")


##find mean of mines and mean of ref
with(mine.ref, tapply(bulk_density, type, mean))
with(mine.ref, tapply(bulk_density, type, sd))

boxplot(bulk_density~type, data = mine.ref)

##make boxplot
##have to change factor names, so I just created a new dataframe to change the names
mf <- mine.ref #new dataframe
levels(mf$type)[levels(mf$type)=="mine"] <- "Mine"
levels(mf$type)[levels(mf$type)=="reference"] <- "Reference"
View(mf)
minevsref_colpalette <- c("#999999", "#009E73")

mine.ref.boxplot <- ggplot(mf, aes(x=site, y=bulk_density)) + 
  xlab("Site") + 
  ylab("Bulk Density " ~ (g/cm^3)) + 
  geom_boxplot() + 
  facet_grid(~general, scale = "free") + 
  geom_point(data=mf, aes(colour = type, shape = type)) + 
  theme_classic() + 
  scale_colour_manual(values=minevsref_colpalette) + 
  scale_shape_manual(values = c(17, 15)) + 
  theme(strip.text = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(.75, .25),
        legend.justification = c("right", "top"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        text = element_text(size = 10),
        legend.text = element_text(size = 8))
mine.ref.boxplot

library(dplyr)
prot1 = distinct(mf, site, general) %>%
  arrange(site, general)
prot1$yloc = max(mf$bulk_density) + 0.1
prot1$label = c("a", "c", "b", "ab", "ab", "ab", "ab", "ab", "c")
fig2_BD <- mine.ref.boxplot + 
  geom_text(data = prot1, aes(y = yloc, label = label), 
            position = position_dodge(width = .75), size = 4)
fig2_BD

ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\Figure2_BulkDensity.jpeg", plot=fig2_BD, height=11, width=9, units=c("cm"), dpi=600)

##Linear models and test for normality
DVpre.DVR_mod <- lm(bulk_density~type, data=DVpre.DVR)
shapiro.test(resid(DVpre.DVR_mod))
par(mfrow = c(1,2))
plot(DVpre.DVR_mod, which = c(1,2))

##Hines Hill
HH.HHR_mod <- lm(bulk_density~type, data=HH.HHR)
shapiro.test(resid(HH.HHR_mod))
plot(HH.HHR_mod,which=c(1,2))

##SQ
SQpre.SR_mod <- lm(bulk_density~type, data=SQpre.SQR)
shapiro.test(resid(SQpre.SR_mod))
plot(SQpre.SR_mod,which=c(1,2))


###start analysis for pre-rip vs reference

##Two-sample t-tests for mines versus reference
t.test(bulk_density~type, data = DVpre.DVR, var.equal = TRUE)

t.test(bulk_density~type, data = HH.HHR, var.equal = TRUE)

t.test(bulk_density~type, data = SQpre.SQR, var.equal = TRUE)


# Test for differences in variance between mine and reference within each site
library(car)
leveneTest(bulk_density~type, data = DVpre.DVR)
leveneTest(bulk_density~type, data = HH.HHR)
#calculating variances for HH since they are significantly different
aggregate(HH.HHR$bulk_density, by=list(HH.HHR$type), FUN=var)
leveneTest(bulk_density~type, data = SQpre.SQR)
leveneTest(bulk_density~site, data = CT.RR.CTR)
#calculating variances for CT/RR since they are significantly different
aggregate(CT.RR.CTR$bulk_density, by=list(CT.RR.CTR$site), FUN=var)

leveneTest(bulk_density~type, data=mine.ref)


##ANOVA for CT.RR.CRR
CTRRCTR_ANOVA <- aov(bulk_density~site, data = CT.RR.CTR)
shapiro.test(resid(CTRRCTR_ANOVA))
plot(CTRRCTR_ANOVA, which = c(1,2))
summary(CTRRCTR_ANOVA)
TukeyHSD(CTRRCTR_ANOVA)
aggregate(CT.RR.CTR$bulk_density,by=list(CT.RR.CTR$site),FUN=mean)


###stop analysis for pre-rip vs reference


###start analysis for pre-rip vs post-rip number

##SQ flat vs slope
View(SQ.SQR)
wilcox.test(bulk_density~flat_slope, data = SQ.SQR)
library(rstatix)
krustal <- paste(SQ.SQR$flat_slope, SQ.SQR$rip)
SQ.SQR <- cbind(SQ.SQR, krustal)
kruskal.test(bulk_density~krustal, data = SQ.SQR)
wilcox_test(bulk_density~krustal, data = SQ.SQR)


##check for outliers
library(outliers)
test <- grubbs.test(DV.SQ.PrePost$bulk_density)
test


##I changed the lme model to lmer to incorporate the nested "site" variable
pre.post.LMER <- lmer(bulk_density~number*site + (1 | sample), data = DV.SQ.PrePost)
##checking assumptions
plot(pre.post.LMER)
qqnorm(resid(pre.post.LMER))
anova(pre.post.LMER)
summary(pre.post.LMER)
emm2=emmeans(pre.post.LMER, ~number*site)
pairs(emm2, adjust = "holm")
df <- data.frame(pairs(emm2, adjust = "holm"))

##test for dover outlier
hist(DV.SQ.PrePost$bulk_density)
boxplot.stats(DV.SQ.PrePost$bulk_density)$out

# Tests for significant differences in variance across the site/rip combinations
leveneTest(bulk_density~number*site,DV.SQ.PrePost)

library("MuMIn")
r.squaredGLMM(pre.post.LMER)

##boxplot of DV and SQ pre.post.rip
##rename the "number" variables
dspp <- DV.SQ.PrePost #new dataframe
levels(dspp$number)[levels(dspp$number)=="0"] <- "Non-rip"
levels(dspp$number)[levels(dspp$number)=="1"] <- "Single-rip"
levels(dspp$number)[levels(dspp$number)=="2"] <- "Cross-rip"
levels(dspp$number)[levels(dspp$number)=="x"] <- "Pre-rip"
View(dspp)

dv_dspp <- subset(dspp, site=="DV")
with(dv_dspp, tapply(bulk_density, number, mean))
with(dv_dspp, tapply(bulk_density, number, sd))

sq_dspp <- subset(dspp, site=="SQ")
with(sq_dspp, tapply(bulk_density, number, mean))
with(sq_dspp, tapply(bulk_density, number, sd))


library(tidyverse)

pre.post.boxplot <- dspp %>% arrange(number) %>% 
  mutate(number = factor(number, levels = c("Pre-rip", "Non-rip", "Single-rip", "Cross-rip"))) %>% 
  ggplot(aes(x=number, y=bulk_density)) + 
  geom_boxplot(outlier.shape = 3) + 
  xlab("Sample Type") + 
  ylab("Bulk Density"~(g/cm^3)) + 
  facet_wrap(~site) + 
  geom_point(aes(x=number, y=bulk_density, shape = number)) + 
  scale_shape_manual(values = c(15, 16, 17, 3)) + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.position = c(.92, .42),
        legend.justification = c("right", "top"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        text = element_text(size = 10),
        legend.text = element_text(size = 8))
pre.post.boxplot

##add post hoc test
library(dplyr)
prot2 = distinct(dspp, site, number) %>%
  arrange(site, number)
prot2$yloc = max(dspp$bulk_density) + 0.1
prot2$label = c("a", "b", "bc", "b", "bc", "cd", "d", "bc")
fig5_BD <- pre.post.boxplot + 
  geom_text(data = prot2, aes(y = yloc, label = label), 
            position = position_dodge(width = .75), size = 4)
  
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\Figure5_BulkDensity_ripnumber.jpeg", plot=fig5_BD, height=11, width=9, units=c("cm"), dpi=600)
