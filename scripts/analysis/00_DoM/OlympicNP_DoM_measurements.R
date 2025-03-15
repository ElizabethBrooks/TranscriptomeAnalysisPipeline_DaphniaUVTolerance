# script to plot Olympic NP pond DoM measurements

# install nationalparkcolors package
#devtools::install_github("katiejolly/nationalparkcolors")

# import libraries
library(ggplot2)
library(nationalparkcolors)

# set working directory
setwd("/Users/bamflappy/PfrenderLab/OLYM22/")

# import pond data
domData <- read.csv("OLYM_06to22_DoM_report_OLYM_overall_noOnlyDry.csv", header=TRUE)

# retrieve the vector of colors associated with Zion
park_colors <- park_palette("Zion")

# plot Olympic NP Pond Absorbance and Attenuation
ggplot(data=domData, 
       aes(x=A440, y=Kd, group=factor(Year))) +
  geom_point(aes(color=factor(Year), shape=factor(Lake))) +
  theme_bw() +
  #geom_text(aes(label=Primenet)) +
  scale_color_manual(values = park_colors) +
  labs(title = "Olympic NP Pond Absorbance and Attenuation", 
       x ="Absorbance (A440)", 
       y = "Attenuation (Kd)") +
  theme(
    plot.title = element_text(color = park_colors[1], size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.x = element_text(color = park_colors[4], size = 14, face = "bold"),
    axis.title.y = element_text(color = park_colors[5], size = 14, face = "bold")
)

# save the last plot
ggsave("OLYM_06to22_Kd_A440.png", bg = "white")

# plot Olympic NP Pond Transparency (% Surface UVR at 10cm Depth)
ggplot(data=domData, 
       aes(x=reorder(Primenet, pUV10cm), y=pUV10cm, group=factor(Year))) +
  geom_point(aes(color=factor(Year), shape=factor(Lake))) +
  theme_bw() +
  #geom_text(aes(label=Primenet)) +
  scale_color_manual(values = park_colors) +
  labs(title = "Olympic NP Pond Transparency", 
       x ="Ponds by Increasing Transparency", 
       y = "Transparency") +
  theme(
    plot.title = element_text(color = park_colors[1], size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.x = element_text(color = park_colors[4], size = 14, face = "bold"),
    axis.title.y = element_text(color = park_colors[5], size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

# save the last plot
ggsave("OLYM_06to22_transparency.png", bg = "white")
