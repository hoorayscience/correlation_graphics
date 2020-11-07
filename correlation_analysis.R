##Code to assumption check and analyze regression##
##Ellie Wenger##
##11/5/2020##

##Load Packages##
library(dplyr)
library(ggpubr)
library(readr)
library(rstatix)
library(tidyverse)
library(PerformanceAnalytics)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(pander)

##Load tidy data###
PAH_Data <- read_csv("PAH_Data.csv")
set.seed(1234)
dplyr::sample_n(PAH_Data, 10)

numerical.vars <- PAH_Data %>% select_if(is.numeric)
head(numerical.vars)

##check normality of numerical data frame
##still need to figure out how to display results with column names
normality <- list()
for(i in seq(dim(numerical.vars)[2])) {
  normality <- append(normality, shapiro.test(numerical.vars[[i]]))
}
normality

##log transform non-normal variables
trans.numerical.vars <- log(numerical.vars[1:16])
colnames(trans.numerical.vars)

#check normality on log transformed variables
log.normality <- list()
for(i in seq(dim(trans.numerical.vars)[2])) {
  log.normality <- append(log.normality, shapiro.test(trans.numerical.vars[[i]]))
}
log.normality

#If some variables are over transformed, replace the column in the log data frame wih a root transformation
trans.numerical.vars$`Benzo(k)fluoranthene` <- numerical.vars$`Benzo(k)fluoranthene`^(1/4)
trans.numerical.vars$`Benz(a)anthracene` <- numerical.vars$`Benz(a)anthracene`^(1/4)

##use correlation chart to check out histograms and pearson correlations
correlation <- trans.numerical.vars[, c(1:16)]
chart.Correlation(correlation, histogram=TRUE, method=c("pearson"), pch=19)

##for loop for correlation tests. does not work yet

for (i in 1:length(trans.numerical.vars)) {
  a <- cor.test(trans.numerical.vars$`AhR BEQ (ng TCDD/g)`, trans.numerical.vars[,i], method=c("pearson"), use = "complete.obs")
  print(paste(colnames(trans.numerical.vars)[i], " est:", a$estimate, " p=value:", a$p.value))
}

##for loop to plot all congeners

for(i in 9:length(PAH_Data)) {                              # ggplot within for-loop
  ggplot(PAH_Data, aes(x = AhR_BEQ_ng_TCDD_g, y = PAH_Data[ , i])) +
    geom_point(size=3, aes(colour=Stratum, shape=Stratum))+
    geom_smooth(method = "lm", colour="black", lwd=0.5) +
    ylim(0,10) +
    ylab ("PCB118") +
    theme_classic() +
    theme(legend.position = "none") + 
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  Sys.sleep(2)
}

#maual ggplot instead of loop ): I got a deadline
Total_PAHs_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Total_PAHs)) +
  geom_point(size=3, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Total PAH") +
  xlab("AhR BEQ (ng TCDD/g)" )+
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12))

Total_PAHs_AhR_graph

#Plot each congener,
#had a lot of trouble plotting just one trendline
#The 'trick' is to specify the grouping variable in geom_point only, and not in the general aesthetics.

Benzo.b.fluoranthene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Benzo.b.fluoranthene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Benzo(b)fluoranthene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Benzo.b.fluoranthene_AhR_graph

Benzo.k.fluoranthene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Benzo.k.fluoranthene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Benzo(k)fluoranthene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Benzo.k.fluoranthene_AhR_graph

Dibenz.a.h.anthracene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Dibenz.a.h.anthracene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Dibenz(a,h)anthracene") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank())

Dibenz.a.h.anthracene_AhR_graph

Indeno.1.2.3.cd.pyrene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Indeno.1.2.3.cd.pyrene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Indeno(1,2,3-cd)pyrene") +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Indeno.1.2.3.cd.pyrene_AhR_graph

Benzo.a.pyrene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Benzo.a.pyrene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Benzo(a)pyrene") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Benzo.a.pyrene_AhR_graph

Benzo.e.pyrene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Benzo.e.pyrene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Benzo(e)pyrene") +
  xlab("AhR BEQ (ng TCDD/g)")+
  theme_classic() + 
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Benzo.e.pyrene_AhR_graph

Benz.a.anthracene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Benz.a.anthracene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  ylab ("Benz(a)anthracene") +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Benz.a.anthracene_AhR_graph

Chrysene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Chrysene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  xlab("AhR BEQ (ng TCDD/g)") +
  ylab ("Chrysene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Chrysene_AhR_graph

Fluoranthene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Fluoranthene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  xlab("AhR BEQ (ng TCDD/g)") +
  ylab ("Fluoranthene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Fluoranthene_AhR_graph

Pyrene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Pyrene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  xlab("AhR BEQ (ng TCDD/g)") +
  ylab ("Pyrene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Pyrene_AhR_graph

Phenanthrene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Phenanthrene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  xlab("AhR BEQ (ng TCDD/g)") +
  ylab ("Phenanthrene") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank())

Phenanthrene_AhR_graph

Anthracene_AhR_graph <- ggplot(PAH_Data, aes(x=AhR_BEQ_ng_TCDD_g, y=Anthracene)) +
  geom_point(size=2.5, aes(colour=Stratum, shape=Stratum)) +
  geom_smooth(method = "lm", colour="black", lwd=0.5) +
  xlab("AhR BEQ (ng TCDD/g)") +
  ylab ("Anthracene") +
  theme_classic() +
  theme(legend.position = "bottom")+
  theme(axis.title.x=element_blank())
        
Anthracene_AhR_graph

#arrange into one graphic
PAH_congener_grid <- grid.arrange(Benzo.b.fluoranthene_AhR_graph, Benzo.k.fluoranthene_AhR_graph, Dibenz.a.h.anthracene_AhR_graph, Indeno.1.2.3.cd.pyrene_AhR_graph, Benzo.a.pyrene_AhR_graph, Benzo.e.pyrene_AhR_graph, Benz.a.anthracene_AhR_graph, Chrysene_AhR_graph, Fluoranthene_AhR_graph, Pyrene_AhR_graph, Phenanthrene_AhR_graph, Anthracene_AhR_graph, ncol = 2)

theme(axis.text.y=element_blank(),
      axis.ticks.y=element_blank())

