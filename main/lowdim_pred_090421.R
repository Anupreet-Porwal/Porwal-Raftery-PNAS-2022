#### Libraries #### ----
library(mlbench)
library(reshape2)
library(ggplot2)
library(MASS)
library(BMA)
library(glmnet)
library(BAS)
library(EMVS)
library(SSLASSO)
library(horseshoe)
library(mombf)
library(networkBMA)
library(ncvreg)
library(BoomSpikeSlab)
library(ModelMetrics)
library(robustHD)
library(xtable)
library(ISLR)
library(parallel)
library(doParallel)
library(plyr)
library(xlsx)
library(XLConnect)
library(openxlsx)
library(SIS)
library(gss)
library(spikeslab)

registerDoParallel(cores=2)

#################### EDIT THIS PART ONLY ######################
#setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")

args = commandArgs(TRUE)

data.num <- as.numeric(args[1])
runtype <- as.character(args[3])


if(runtype=='BMS'){
  source("src/bms_prediction_study_function_083121.R")
}else if(runtype=='BMA'){
  source("src/prediction_study_function_011321.R")
}

datanames <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily",
               "BS-hourly","SML")


dataname <- datanames[data.num]


resultsFile = paste(as.character(args[2]), paste0(dataname,
                                                  "_","prediction","_",runtype),sep = "")


if(dataname=="College"){
  data("College")
  prelim.data=College

  col.num <- ncol(prelim.data)
  prelim.data=prelim.data[ , c(1,3:(col.num),2)]
  # removed enroll and accept due to causual issues
  prelim.data=subset(prelim.data, select=-c(Enroll,Accept))
  prelim.data$Apps=log(prelim.data$Apps)
  prelim.data$F.Undergrad=log(prelim.data$F.Undergrad)
  prelim.data$P.Undergrad=log(prelim.data$P.Undergrad)

}else if(dataname=="BC-Tmax"){
  prelim.data <- read.table("data/BiasCorrection/Bias_correction_ucl.csv",header = TRUE,sep = ",")

  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , -c(1:2,colnum)]
  prelim.data <- prelim.data[complete.cases(prelim.data),]
  prelim.data[ , ncol(prelim.data)] <- sqrt(prelim.data[ , ncol(prelim.data)] )


}else if(dataname=="BC-Tmin"){
  prelim.data <- read.table("data/BiasCorrection/Bias_correction_ucl.csv",header = TRUE,sep = ",")

  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , -c(1:2,colnum-1)]
  prelim.data <- prelim.data[complete.cases(prelim.data),]
  prelim.data[ , ncol(prelim.data)] <- sqrt(prelim.data[ , ncol(prelim.data)] )

}else if(dataname=="BS-daily"){
  prelim.data=read.table("data/bikesharing/day.csv",header = TRUE,sep = ",")

  # remove date and record index since these variables are just indexing variables
  # remove casual and registered since they define the outcome variable of interest
  prelim.data=prelim.data[ ,-c(1,2,14,15)]

  # Changing normalised versions to actual numbers based on the read me file
  prelim.data$temp <- prelim.data$temp*41
  prelim.data$atemp <- prelim.data$atemp*50
  prelim.data$windspeed <- prelim.data$windspeed*67
  prelim.data$hum <- prelim.data$hum*100

  prelim.data$yr <- factor(format(prelim.data$yr, format="%A"),
                           levels = c("0", "1") , labels = c("2011","2012"))
  prelim.data$weathersit <- factor(format(prelim.data$weathersit, format="%A"),
                                   levels = c("1", "2","3") ,
                                   labels = c("Good","Moderate","Bad"))
  prelim.data$holiday <- factor(format(prelim.data$holiday, format="%A"),
                                levels = c("0", "1") , labels = c("NotHoliDay","Holiday"))
  prelim.data$season <- factor(format(prelim.data$season, format="%A"),
                               levels = c("1", "2","3","4") , labels = c("Spring","Summer","Fall","Winter"))
  prelim.data$mnth <- factor(prelim.data$mnth)
  prelim.data$mnth <- relevel(prelim.data$mnth,ref=6)
  prelim.data$weekday <- factor(prelim.data$weekday)
  prelim.data$weekday <- relevel(prelim.data$weekday,ref=3)
  prelim.data=subset(prelim.data,select = -c(workingday))
  prelim.data$cnt=(prelim.data$cnt)^(1/2)

}else if(dataname=="BS-hourly"){

  prelim.data=read.table("data/bikesharing/hour.csv",header = TRUE,sep = ",")

  # remove date and record index since these variables are just indexing variables
  # remove casual and registered since they define the outcome variable of interest
  prelim.data=prelim.data[ ,-c(1,2,15,16)]

  # Changing normalised versions to actual numbers based on the read me file
  prelim.data$temp <- prelim.data$temp*41
  prelim.data$atemp <- prelim.data$atemp*50
  prelim.data$windspeed <- prelim.data$windspeed*67
  prelim.data$hum <- prelim.data$hum*100

  prelim.data$yr <- factor(format(prelim.data$yr, format="%A"),
                           levels = c("0", "1") , labels = c("2011","2012"))
  prelim.data$weathersit <- factor(format(prelim.data$weathersit, format="%A"),
                                   levels = c("1", "2","3","4") ,
                                   labels = c("Good","Moderate","Bad","Bad"))
  prelim.data$holiday <- factor(format(prelim.data$holiday, format="%A"),
                                levels = c("0", "1") , labels = c("NotHoliDay","Holiday"))
  prelim.data$season <- factor(format(prelim.data$season, format="%A"),
                               levels = c("1", "2","3","4") , labels = c("Spring","Summer","Fall","Winter"))
  prelim.data$hr[prelim.data$hr %in% 0:5] <- "LateNight"
  prelim.data$hr[prelim.data$hr %in% 6:8] <- "EarlyMorning"
  prelim.data$hr[prelim.data$hr %in% 9:15] <- "Morning"
  prelim.data$hr[prelim.data$hr %in% 16:18] <- "Evening"
  prelim.data$hr[prelim.data$hr %in% 19:23] <- "Night"

  prelim.data$hr <- as.factor(prelim.data$hr)

  prelim.data$hr <- relevel(prelim.data$hr,ref="Morning")

  # Cube root tranformation of dependent variable, season: ref=summer
  # month: ref=may, weekday:ref=wednesday, hr:ref=Morning

  prelim.data$season <- relevel(prelim.data$season,ref=2)
  prelim.data$mnth <- factor(prelim.data$mnth)
  prelim.data$mnth <- relevel(prelim.data$mnth,ref=5)
  prelim.data$weekday <- factor(prelim.data$weekday)
  prelim.data$weekday <- relevel(prelim.data$weekday,ref=3)
  prelim.data=subset(prelim.data,select = -c(workingday))
  prelim.data$cnt=(prelim.data$cnt)^(1/3)

}else if(dataname=="SML"){

  prelim.data <- read.table("data/SML2010/NEW-DATA-2.T15.txt",header = TRUE,sep = " ")

  prelim.data <- prelim.data[ , -c(1:2,19:21)]
  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , c(3:colnum,1)]
  prelim.data$X24.Day_Of_Week <- round(prelim.data$X24.Day_Of_Week,digits = 0)
  prelim.data$X24.Day_Of_Week <- factor(prelim.data$X24.Day_Of_Week)

}

# Standardizing all the x's except the factor variables
for (i in 1:(ncol(prelim.data)-1)){
  if(is.factor(prelim.data[ ,i])==0){
    prelim.data[ ,i]=(prelim.data[ ,i]-mean(prelim.data[ ,i]))/sd(prelim.data[ ,i])
  }
}

# centering the Y. Assuming Y is always the last column of prelim.data
#prelim.data[ , ncol(prelim.data)]=prelim.data[ , ncol(prelim.data)]-mean(prelim.data[ , ncol(prelim.data)])
p <- ncol(prelim.data)-1
colnames(prelim.data)[1:p] <- paste("V",1:p,sep="")
colnames(prelim.data)[ncol(prelim.data)]="Y"

#### NO EDITS BELOW THIS LINE (except RESULTS) ####
datamat <- prelim.data

Xmat=datamat[ ,-ncol(datamat)]
Y=datamat[ ,ncol(datamat)]

if(runtype=='BMS'){
  res1 <- prediction_study_bms(Xmat = Xmat,Y=Y,bootn = 100,split = 0.75)
}else if(runtype=='BMA'){
  res1 <- prediction_study(Xmat = Xmat,Y=Y,bootn = 100,split = 0.75)
}

write.xlsx(res1,file=paste(resultsFile,"xlsx",sep="."))
save(res1, file = paste(resultsFile,"rda",sep = "."))


