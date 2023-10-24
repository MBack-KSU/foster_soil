##Author: Mike Back
##Organization: Kent State University
##Project: FoSTER Project, 2023 manuscript
##Script purpose: Ksat stats and figures
##Date updated:7/26/2023


##Import data
#call "Ksat.csv"
Ksat <- read.csv("C:\\Users\\mikeb\\Documents\\github\\foster_soil\\data\\ksat.csv", stringsAsFactors = T)
View(Ksat)

library("nlme")
library("lme4")
library("lmerTest")
library("emmeans")
library("ggplot2")
library("MuMIn")
library("scales")

###Create subsets

##Dover
DV.DVR <- subset(Ksat, site=="DV"|site=="DVR")
Dpre <- subset(DV.DVR, pre_post=="pre")
Dpost <- subset(DV.DVR, pre_post=="post")
DR <- subset(Ksat, site=="DVR")
Dpre.DR <- subset(DV.DVR, pre_post=="pre"|site=="DVR")
Dpost.DR <- subset(DV.DVR, pre_post=="post"|site=="DVR")
D.pre.post <- subset(Ksat, site=="DV")

##Snowville
SQ.SQR <- subset(Ksat, site=="SQ"|site=="SQR")
SQpre <- subset(SQ.SQR, pre_post=="pre")
SQpost <- subset(SQ.SQR, pre_post=="post")
SQR <- subset(Ksat, site=="SQR")
SQpre.SQR <- subset(SQ.SQR, pre_post=="pre"|site=="SQR")
SQpost.SQR <- subset(SQ.SQR, pre_post=="post"|site=="SQR")
SQ.pre.post <- subset(Ksat, site=="SQ")

##Dover and SQ
DV.SQ.PrePost <- subset(Ksat, site=="SQ"|site=="DV")

##Hines Hill
HH.HR <- subset(Ksat, site=="HH"|site=="HHR")
HH <- subset(Ksat, site=="HH")
HR <- subset(Ksat, site=="HHR")

##Cleveland Trust, Rockside Road, Cleaveland Trust and Rockside Road Reference
CT.RR.CRR <- subset(Ksat, site=="CT"|site=="RR"|site=="CTR")
CT.CRR <- subset(Ksat, site=="CT"|site=="CTR")
RR.CRR <- subset(Ksat, site=="RR"|site=="CTR")
CT <- subset(Ksat, site=="CT")
RR <- subset(Ksat, site=="RR")
CRR <- subset(Ksat, site=="CTR")

##subset for all mine vs ref
mine.ref <- subset(Ksat, pre_post=="pre"|pre_post=="no")
summary(mine.ref)


##Mixed effects model for all mines vs their respective reference sites, needed log transformation
mine.ref.LMER <- lmer(log(Ksat.ms)~ref_mine + (1 | general), data = mine.ref)
plot(mine.ref.LMER)
qqnorm(resid(mine.ref.LMER))
anova(mine.ref.LMER)
summary(mine.ref.LMER)
emm1=emmeans(mine.ref.LMER, ~ref_mine)
pairs(emm1, adjust = "holm")
r.squaredGLMM(mine.ref.LMER)
##test model without general as a random effect
mine.ref.LM <- lm(log(Ksat.ms)~ref_mine, data = mine.ref)
plot(mine.ref.LM)
qqnorm(resid(mine.ref.LM))
anova(mine.ref.LM)
summary(mine.ref.LM)
##likelihood ratio test to test significance of random effect
anova(mine.ref.LMER, mine.ref.LM)

##linear model to determine significant differences in Ksat between sites
mine.ref.sites <- lm(log(Ksat.ms)~site, data=mine.ref)
anova(mine.ref.sites)
emm3=emmeans(mine.ref.sites, ~site)
pairs(emm3, adjust = "holm")#post hoc test used for significance letters on figure

##find mean for all mine sites combined and references combined
with(mine.ref, tapply(Ksat.ms, ref_mine, mean))
with(mine.ref, tapply(Ksat.ms, ref_mine, sd))

boxplot(log(Ksat.ms)~ref_mine, data = mine.ref)

##FIGURE 2B
#code for building Ksat part of fig 2, combine Ksat and bulk density in powerpoint

#Create new dataframe to rename "mine" and "ref" factors
mf <- mine.ref #new dataframe
levels(mf$ref_mine)[levels(mf$ref_mine)=="mine"] <- "Mine"
levels(mf$ref_mine)[levels(mf$ref_mine)=="ref"] <- "Reference"

#color blind friendly palette
minevsref_colpalette <- c("#999999", "#009E73")

mine.ref.boxplot <- ggplot(mf, aes(x=site, y=Ksat.ms)) + 
  labs(x = "Site", y = "Ksat" ~ (m/s)) +
  geom_boxplot() + 
  facet_grid(~general, scale = "free") + 
  geom_point(data=mf, aes(colour = ref_mine, shape = ref_mine)) + 
  theme_classic() + 
  scale_colour_manual(values=minevsref_colpalette) + 
  scale_shape_manual(values = c(17, 15)) + 
  theme(strip.text = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(.95, .35),
        legend.justification = c("right", "top"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        text = element_text(size = 10),
        legend.text = element_text(size = 8)) + 
  scale_y_log10()
mine.ref.boxplot

#code to add significance letters
library(dplyr)
prot1 = distinct(mf, site, general) %>%
  arrange(site, general)
prot1$yloc = max(mf$Ksat.ms) + 0.0001
prot1$label = c("a", "bc", "bc", "bc", "ab", "bc", "bc", "c", "bc")#RR is not correct (third letter), it should be an "a", fix in powerpoint
fig2_ksat <- mine.ref.boxplot + 
  geom_text(data = prot1, aes(y = yloc, label = label), 
            position = position_dodge(width = .75), size = 4)

fig2_ksat

ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\Figure2_Ksat.jpeg", plot=fig2_ksat, height=11, width=9, units=c("cm"), dpi=600)


###Linear models and tests for normality

##Lm SQ vs SQR
SQpre.SQR_mod <- lm(Ksat.ms~ref_mine, data = SQpre.SQR)
plot(SQpre.SQR_mod, which = c(1,2))

##Transform Snowville Ksat variable with a log function
SQpre.SQR$logKsat <- log(SQpre.SQR$Ksat.ms)

##Lm SQ vs SQR, using log(Ksat)
SQpre.SQR_logmod <- lm(logKsat~ref_mine, data = SQpre.SQR)
shapiro.test(resid(SQpre.SQR_logmod))
plot(SQpre.SQR_logmod, which = c(1,2))

##Linear model for Hines Hill versus Reference
HH.HR_mod <- lm(Ksat.ms~ref_mine, data = HH.HR)
shapiro.test(resid(HH.HR_mod))
plot(HH.HR_mod, which = c(1,2))

##Linear model for Dover pre-rip versus Dover Reference
Dpre.DR_mod <- lm(Ksat.ms~ref_mine, data = Dpre.DR)
shapiro.test(resid(Dpre.DR_mod))
plot(Dpre.DR_mod, which = c(1,2))


###Start ref vs mine analysis

##Two sample t tests

t.test(Ksat.ms~ref_mine, data = HH.HR, var.equal = TRUE)

t.test(Ksat.ms~ref_mine, data = Dpre.DR, var.equal = TRUE)

t.test(log(Ksat.ms)~ref_mine, data = SQpre.SQR, var.equal = TRUE)


# Levene's test for unequal variances
library(car)
leveneTest(Ksat.ms~ref_mine, data = HH.HR)
leveneTest(Ksat.ms~ref_mine, data = Dpre.DR)
leveneTest(Ksat.ms~ref_mine, data = SQpre.SQR)

##ANOVA model for CT, RR, and CRR
CT.RR.CRR_ANOVA <- aov(log(Ksat.ms)~site, data = CT.RR.CRR)
shapiro.test(resid(CT.RR.CRR_ANOVA))
plot(CT.RR.CRR_ANOVA, which = c(1,2))
summary(CT.RR.CRR_ANOVA)
TukeyHSD(CT.RR.CRR_ANOVA)


##log trans used above seems to fit better than the asinsqrt she used here below  -mike
##Transform CT.RR.CRR data using asinsqrt
CT.RR.CRR$asinKsat.ms <- asin(sqrt(CT.RR.CRR$Ksat.ms))

##ANOVA model for CT.RR.CRR using asinsqrt transformation
CTRRCRRasin_ANOVA <- aov(asinKsat.ms~site, data = CT.RR.CRR)
shapiro.test(resid(CTRRCRRasin_ANOVA))
plot(CTRRCRRasin_ANOVA, which = c(1,2))
summary(CTRRCRRasin_ANOVA)
TukeyHSD(CTRRCRRasin_ANOVA)

###Stop ref vs mine analysis


### Start analysis for Flat vs slope and rip orientation
View(SQ.pre.post)

##Lm SQ Flat vs Slope
SQ.flat.slope.mod <- lm(Ksat.ms~F_S, data = SQ.pre.post)
shapiro.test(resid(SQ.flat.slope.mod))
plot(SQ.flat.slope.mod, which = c(1,2))

##Lm SQ FLat vs Slope, using log trans
SQ.pre.post$logKsat <- log(SQ.pre.post$Ksat.ms)
SQ.flat.slope.logmod <- lm(logKsat~F_S, data = SQ.pre.post)
shapiro.test(resid(SQ.flat.slope.logmod))
plot(SQ.flat.slope.logmod, which = c(1,2))

##t test with log mod for flat vs slope
t.test(logKsat~F_S, data = SQ.pre.post, var.equal = TRUE)

##SQ Flat vs Slope, using non-parametric ranks (Mann-Whitney-U test)

wilcox.test(Ksat.ms~F_S, data = SQ.pre.post)
library(rstatix)
krustal <- paste(SQ.pre.post$F_S, SQ.pre.post$pre_post)
SQ.pre.post <- cbind(SQ.pre.post, krustal)
kruskal.test(Ksat.ms~krustal, data = SQ.pre.post)
wilcox_test(Ksat.ms~krustal, data = SQ.pre.post)

##Lm SQ orientation
SQ.orientation.mod <- lm(log(Ksat.ms)~orientation, data = SQ.pre.post)
shapiro.test(resid(SQ.orientation.mod))
plot(SQ.orientation.mod, which = c(1,2))

##t test with log mod for orientation
t.test(log(Ksat.ms)~orientation, data = SQ.pre.post, var.equal = TRUE)

##SQ orientation, using using non-parametric ranks (Mann-Whitney-U test)
wilcox.test(Ksat.ms~orientation, data = SQ.pre.post)

krustal2 <- paste(SQ.pre.post$orientation, SQ.pre.post$pre_post)
SQ.pre.post <- cbind(SQ.pre.post, krustal2)
kruskal.test(orientation~krustal2, data = SQ.pre.post)
wilcox_test(orientation~krustal2, data = SQ.pre.post)

###Stop analysis for Flat vs slope and rip orientation


###Start pre vs post analysis

##I'm just using the same code as we did for bulk density
pre.post.LMER <- lmer(log(Ksat.ms)~rip_number*site + (1 | sample), data = DV.SQ.PrePost)
##checking assumptions
plot(pre.post.LMER)
qqnorm(resid(pre.post.LMER))
qqline(resid(pre.post.LMER),col="red")
plot(fitted(pre.post.LMER),resid(pre.post.LMER))
anova(pre.post.LMER)
summary(pre.post.LMER)
emm2=emmeans(pre.post.LMER, ~rip_number*site)
pairs(emm2, adjust = "holm")

# Levene's test for unequal variances
leveneTest(log(Ksat.ms)~rip_number*site, data = DV.SQ.PrePost)
aggregate(log(DV.SQ.PrePost$Ksat.ms), by=list(paste(DV.SQ.PrePost$rip_number,DV.SQ.PrePost$site)), FUN=var)

##boxplot but have to rename "rip_number" variables
dsppksat <- DV.SQ.PrePost #new dataframe
levels(dsppksat$rip_number)[levels(dsppksat$rip_number)=="0"] <- "Non-rip"
levels(dsppksat$rip_number)[levels(dsppksat$rip_number)=="1"] <- "Single-rip"
levels(dsppksat$rip_number)[levels(dsppksat$rip_number)=="2"] <- "Cross-rip"
levels(dsppksat$rip_number)[levels(dsppksat$rip_number)=="X"] <- "Pre-rip"

pre.post.boxplot <- ggplot(dsppksat, aes(x=rip_number, y=Ksat.ms)) + 
  geom_boxplot() + 
  xlab("Sample Type") + 
  ylab("Ksat" ~ (m/s)) + 
  facet_wrap(~site) + 
  scale_y_log10() + 
  geom_point(data=dsppksat, aes(x=rip_number, y=Ksat.ms, shape = rip_number)) + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.position = c(.35, .405),
        legend.justification = c("right", "top"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 10),
        legend.text = element_text(size = 8))
pre.post.boxplot

ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\Figure5_Ksat_ripnumber.jpeg", plot=pre.post.boxplot, height=11, width=9, units=c("cm"), dpi=600)
