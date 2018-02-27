## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----loading-libraries-and-data------------------------------------------
library(lori)
library(FactoMineR)
library(logmult)
death <- read.table("http://factominer.free.fr/book/death.csv", header = TRUE,
                    sep = ";", row.names = 1)
colnames(death) <- c("0-1", "1-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74",
                    "75-84", "85-94", "95 and more")
## Extract first 65 rows correponding to year 2006
death <- death[66:130,]
rownames(death) <- sapply(rownames(death), function(s) gsub("^.*?_","",s))
n <- nrow(death)
p <- ncol(death)

## ----apply-lori, echo=FALSE----------------------------------------------

## Run the algo 
death_fit_lori <- lori(death)

## Extract model parameters
mu_lori <- death_fit_lori$mu
alpha_lori <- death_fit_lori$alpha
beta_lori <- death_fit_lori$beta
theta_lori <- death_fit_lori$Theta
colnames(theta_lori) <- colnames(death)
rownames(theta_lori) <- rownames(death)

## ----display-plots-------------------------------------------------------
plot_interaction(theta_lori, title = "Display plot of LORI model", selectRow = "contrib 10")

