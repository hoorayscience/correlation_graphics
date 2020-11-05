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

correlation <- trans.numerical.vars[, c(1:16)]
chart.Correlation(correlation, histogram=TRUE, pch=19)




