############################################################################################

############################################################################################
# Setting up a working directory:
setwd('~/Desktop')


# install.packages("devtools")

library(ggplot2)
library(stargazer)
library(sjPlot)
library(MASS)
library(plyr)

#############################################

#############################################
# Load Dataset (update this path to the csv file's location)
mydata = read.csv("///CerebellarQMRI_data.csv")

mydata$Age_Group <- as.factor(mydata$Age_Group)

# Define age groups:
Children <- mydata[which(mydata$Age<17),]
Adults <- mydata[which(mydata$Age>17),]


# Subtract from each lobule the childhood mean (adult distribution thus varies relative to this point)
# Children lobule means: 
childhood_meanT1_I_IV <- mean(Children$meanT1_I_IV, na.rm=T)
childhood_meanT1_Crus1 <- mean(Children$meanT1_Crus1, na.rm=T)

childhood_meanTV_I_IV <- mean(Children$meanTV_I_IV, na.rm=T)
childhood_meanTV_Crus1 <- mean(Children$meanTV_Crus1, na.rm=T)


# Subtract childhood means from all data
mydata$childhood_centered_meanT1_I_IV <- (mydata$meanT1_I_IV - childhood_meanT1_I_IV)
mydata$childhood_centered_meanT1_Crus1 <- (mydata$meanT1_Crus1 - childhood_meanT1_Crus1)

# Subtract childhood means from all data
mydata$childhood_centered_meanTV_I_IV <- (mydata$meanTV_I_IV - childhood_meanTV_I_IV)
mydata$childhood_centered_meanTV_Crus1 <- (mydata$meanTV_Crus1 - childhood_meanTV_Crus1)

mydata2 <- subset(mydata, select= -c(meanT1_I_IV, meanT1_Crus1, meanTV_I_IV, meanTV_Crus1))

#############################################

#############################################
# Flipped x-axis, and split violin plot:
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#############################################

#############################################

# Make Plot For T1 Relaxation:
library(tidyverse)

#mydata.gatheredT1 <- mydata2 %>%
#	as_data_frame() %>%
#	gather(key = "variable", value="value", -SubjID, -Age, -Age_Group, -Index, -childhood_centered_meanTV_I_IV, -childhood_centered_meanTV_Crus1)


#Relevel Labels
#mydata.gatheredT1$variable <- factor(mydata.gatheredT1$variable, levels=c("childhood_centered_meanT1_I_IV",  "childhood_centered_meanT1_Crus1"))

#Plot
#ggplot(na.omit(mydata.gatheredT1), aes(x= variable, y= value, fill=Age_Group)) + theme_minimal() + labs(title="", x="", y="\nT1 Relaxation [sec]") + scale_fill_discrete(name="Age_Group", labels=c("Children", "Adults")) + theme(axis.text = element_text(size = 11), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(size = 11), legend.position = "right", legend.title = element_blank()) + scale_x_discrete(labels=c("Lobule I-IV", "Crus I")) + geom_split_violin() + geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0, show.legend=F) + theme(legend.key.size=unit(0.5, "cm"))


# Make Plot For MTV:
mydata.gatheredTV <- mydata2 %>%
	as_data_frame() %>%
	gather(key = "variable", value="value", -SubjID, -Age, -Age_Group, -Index, -childhood_centered_meanT1_I_IV,  -childhood_centered_meanT1_Crus1)

#Relevel Labels
mydata.gatheredTV$variable <- factor(mydata.gatheredTV$variable, levels=c("childhood_centered_meanTV_I_IV",  "childhood_centered_meanTV_Crus1"))


#Plot
ggplot(na.omit(mydata.gatheredTV), aes(x= variable, y= value, fill=Age_Group)) + theme_minimal() + labs(title="", x="", y="\nMacromolecular Tissue Volume\n[milliliters, normalized to children]") + scale_fill_discrete(name="Age Group", labels=c("Children", "Adults")) + theme(axis.text = element_text(size = 11), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), axis.title = element_text(size = 11), legend.position = "right", legend.title = element_blank()) + scale_x_discrete(labels=c("Lobule I-IV",  "Crus I")) + geom_split_violin() + geom_boxplot(width = 0.25, notch = FALSE, notchwidth = .4, outlier.shape = NA, coef=0, show.legend=F) + theme(legend.key.size=unit(0.5, "cm"))

