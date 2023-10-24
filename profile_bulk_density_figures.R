##Author: Mike Back
##Organization: Kent State University
##Project: FoSTER Project, 2023 manuscript
##Script purpose: Building profile bulk density figures
##Date updated:7/26/2023


##load data
profiles <- read.csv("C:\\Users\\mikeb\\Documents\\github\\foster_soil\\data\\profile_bulk_density.csv", stringsAsFactors = T)
summary(profiles)

##instal packages
library(tidyverse)
library(scales)
library(ggpubr)

View(SQ.flat.pre_rip)
###time to subset
#these change everytime the data file changes...so dont change the csv...

##SQ ref vs pre
SQ.ref <- profiles[c(65:76), c(1:5)]
SQ.slope.pre_rip <- profiles[c(85:93), c(1:5)]
SQ.flat.pre_rip <- profiles[c(77:84), c(1:5)]

##HH ref vs pre
HH.ref <- profiles[c(94:102), c(1:5)]
HH.pre <- profiles[c(103:114), c(1:5)]

##Dover
DV.rip <- profiles[c(115:126), c(1:5)]
DV.non_rip <- profiles[c(127:138), c(1:5)]
View(SQ.slopeEW.non_rip)
##flat
SQ.flat.non_rip <- profiles[c(1:7), c(1:5)]
SQ.flat.rip <- profiles[c(14:20), c(1:5)]
SQ.flat.pre_rip <- profiles[c(77:84), c(1:5)]

##slopeN-S
SQ.slope.non_rip <- profiles[c(21:30), c(1:5)]
SQ.slope.rip <- profiles[c(31:40), c(1:5)]
SQ.slope.pre_rip <- profiles[c(85:93), c(1:5)]

##slope_E-W
SQ.slopeEW.non_rip <- profiles[c(41:52), c(1:5)]
SQ.slopeEW.rip <- profiles[c(53:64), c(1:5)]
SQ.slopeEW.pre_rip <- profiles[c(85:93), c(1:5)]

##palette colors for SQ. Colorblind friendly.
SQ.pre.vs.ref.palette <- c("#999999", "#009E73")
HH.palette <- c("#999999", "#009E73")
DV.palette <- c("#E69F00", "#56B4E9")
SQ.post.rip.Palette <- c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##SQ ref vs pre
a <- ggplot(SQ.ref, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Reference"), shape = 15) + 
  geom_line(aes(color = "Reference")) + 
  geom_line(data = SQ.flat.pre_rip, aes(color = "Pre-rip")) +  
  geom_point(data = SQ.flat.pre_rip, aes(color = "Pre-rip"), shape = 17) + 
  geom_line(data = SQ.slope.pre_rip, aes(color = "Pre-rip")) + 
  geom_point(data = SQ.slope.pre_rip, aes(color = "Pre-rip"), shape = 17) + 
  ylim(0.9,2) + 
  geom_segment(aes(y = 1.5615660, x = 30, yend = 1.4551698, xend = 40), color = "#999999") + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab("Depth" ~(cm)) + 
  ylab(NULL) + 
  ggtitle("SQ") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") + 
  scale_colour_manual(values=SQ.pre.vs.ref.palette)+ 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))

##HH ref vs pre
b <- ggplot(HH.ref, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Reference"), shape = 15) + 
  geom_line(aes(color = "Reference")) + 
  geom_line(data = HH.pre, aes(color = "Pre-rip")) +  
  geom_point(data = HH.pre, aes(color = "Pre-rip"), shape = 17) + 
  ylim(0.9,2) + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab(NULL) + 
  ylab(NULL) + 
  ggtitle("HH") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") +
  scale_colour_manual(values=HH.palette)+ 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))


##Dover
c <- ggplot(DV.rip, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Rip")) + 
  geom_line(aes(color = "Rip")) + 
  geom_line(data = DV.non_rip, aes(color = "Non-rip")) +  
  geom_point(data = DV.non_rip, aes(color = "Non-rip")) + 
  annotate("segment", x = 30, xend = 35, y = 1.341103, yend = 2, color = "#E69F00", linetype = "dashed") + 
  annotate("point", x = 35, y = 2, color = "#E69F00") + 
  annotate("segment", x = 35, xend = 50, y = 2, yend = 2, color = "#E69F00", linetype = "dashed") + 
  annotate("point", x = 40, y = 2, color = "#E69F00") + 
  annotate("point", x = 45, y = 2, color = "#E69F00") + 
  annotate("point", x = 50, y = 2, color = "#E69F00") + 
  annotate("segment", x = 50, xend = 55, y = 2, yend = 1.617036, color = "#E69F00", linetype = "dashed") + 
  ylim(0.9,2) + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab(NULL) + 
  ylab(NULL) + 
  ggtitle("DV") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") +
  scale_colour_manual(values=DV.palette)+ 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))


##flat
d <- ggplot(SQ.flat.non_rip, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Non-rip")) + 
  geom_line(aes(color = "Non-rip")) + 
  geom_line(data = SQ.flat.rip, aes(color = "Rip")) + 
  geom_point(data = SQ.flat.rip, aes(color = "Rip")) + 
  geom_line(data = SQ.flat.pre_rip, aes(color = "Pre-rip")) + 
  geom_point(data = SQ.flat.pre_rip, aes(color = "Pre-rip"), shape = 17) + 
  ylim(0.9,2) + 
  geom_segment(aes(y = 1.5615660, x = 30, yend = 1.4551698, xend = 40), color = "#999999") + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab("Depth" ~(cm)) + 
  ylab("Bulk Density" ~(g/cm^3)) + 
  ggtitle("SQ Flat") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") + 
  scale_colour_manual(values=SQ.post.rip.Palette)+ 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))

##slopeN-S
e <- ggplot(SQ.slope.non_rip, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Non-rip")) + 
  geom_line(aes(color = "Non-rip")) + 
  geom_line(data = SQ.slope.rip, aes(color = "Rip")) +  
  geom_point(data = SQ.slope.rip, aes(color = "Rip")) + 
  geom_line(data = SQ.slope.pre_rip, aes(color = "Pre-rip")) + 
  geom_point(data = SQ.slope.pre_rip, aes(color = "Pre-rip"), shape = 17) + 
  ylim(0.9,2) + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab(NULL) + 
  ylab("Bulk Density" ~(g/cm^3)) + 
  ggtitle("SQ Slope (N-S)") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") + 
  scale_colour_manual(values=SQ.post.rip.Palette) + 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))

##slopeE-W
f <- ggplot(SQ.slopeEW.non_rip, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Non-rip")) + 
  geom_line(aes(color = "Non-rip")) + 
  annotate("segment", x = 20, xend = 25, y = 1.213384, yend = 2, color = "#E69F00", linetype = "dashed") + 
  annotate("point", x = 25, y = 2, color = "#E69F00") + 
  annotate("segment", x = 25, xend = 45, y = 2, yend = 2, color = "#E69F00", linetype = "dashed") + 
  annotate("point", x = 30, y = 2, color = "#E69F00") + 
  annotate("point", x = 35, y = 2, color = "#E69F00") + 
  annotate("point", x = 40, y = 2, color = "#E69F00") + 
  annotate("point", x = 45, y = 2, color = "#E69F00") + 
  annotate("segment", x = 45, xend = 50, y = 2, yend = 1.450466, color = "#E69F00", linetype = "dashed") + 
  geom_line(data = SQ.slopeEW.rip, aes(color = "Rip")) +  
  geom_point(data = SQ.slopeEW.rip, aes(color = "Rip")) + 
  geom_line(data = SQ.slopeEW.pre_rip, aes(color = "Pre-rip")) + 
  geom_point(data = SQ.slopeEW.pre_rip, aes(color = "Pre-rip"), shape = 17) + 
  ylim(0.9,2) + 
  coord_flip() + 
  scale_x_reverse(breaks = seq(0,60,10)) + 
  xlab(NULL) + 
  ylab("Bulk Density" ~(g/cm^3)) + 
  ggtitle("SQ Slope (E-W)") + 
  geom_hline(yintercept = 1.55, linetype = "dashed", color = "light gray") + 
  scale_colour_manual(values=SQ.post.rip.Palette)+ 
  theme_pubr() +
  theme(legend.position = "none", plot.title = element_text(size = 10), 
        text = element_text(size = 10))

##make a random plot just to get the legend...then cut and paste legend onto new graph
legend_Pal <- c("#E69F00", "#999999", "#009E73", "#56B4E9")

legend_only_please_thanks <- ggplot(HH.ref, aes(x = Depth, y = Bulk_Density)) + 
  geom_point(aes(color = "Reference")) + 
  geom_line(aes(color = "Reference")) + 
  geom_line(data = HH.pre, aes(color = "Pre-rip")) +  
  geom_point(data = HH.pre, aes(color = "Pre-rip")) + 
  geom_point(data = SQ.slope.non_rip, aes(color = "Non-rip")) + 
  geom_line(data = SQ.slope.non_rip, aes(color = "Non-rip")) + 
  geom_line(data = SQ.slope.rip, aes(color = "Rip")) +  
  geom_point(data = SQ.slope.rip, aes(color = "Rip")) + 
  theme_pubr() + 
  scale_colour_manual(values=legend_Pal) + 
  scale_shape_manual(values = c("Reference" = 15, "Pre-rip" = 17)) + 
  theme(legend.position = "right")
legend_only_please_thanks
##add the different shapes in powerpoint cause i can't figure out how to do it in here :(


##combine the 3
fig <- ggarrange(a, b, c, d, e, f, 
                 labels = c("A", "B", "C", "D", "E", "F"), 
                 nrow = 2, ncol = 3,
                 common.legend = T) #crop the legend out in powerpoint
fig
ggsave("G:\\.shortcut-targets-by-id\\0B2rU5M8IUMRwQUd0OE0wT05CR0U\\CVNP data\\Manuscript-BulkDensity-Ksat\\Figure4_profile_BulkDensity.jpeg", plot=fig, height=24, width=19, units=c("cm"), dpi=600)

