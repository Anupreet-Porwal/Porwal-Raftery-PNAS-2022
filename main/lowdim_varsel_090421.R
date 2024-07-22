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
library(xtable)
library(leaps)
library(MAVE)
library(plyr)
library(robustHD)
library(caret)
library(ISLR)
library(parallel)
library(doParallel)
library(car)
library(data.table)
library(openxlsx)
library(SIS)
library(spikeslab)
library(gss)

registerDoParallel(cores=2)

#####
#setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")
#source("src/est_varselection_function_011321.R")

args = commandArgs(TRUE)

data.num <- as.numeric(args[1])
runtype <- as.character(args[3])


if(runtype=='BMS'){
  source("src/bms_est_varselection_function_083121.R")
}else if(runtype=='BMA'){
  source("src/est_varselection_function_011321.R")
}



datanames <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily",
               "BS-hourly","SML")


dataname <- datanames[data.num]


resultsFile = paste(as.character(args[2]), paste0(dataname,
                                                  "_", runtype),sep = "")

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
y.val.0 <- mean(prelim.data[ , ncol(prelim.data)])
prelim.data[ , ncol(prelim.data)]=prelim.data[ , ncol(prelim.data)]-mean(prelim.data[ , ncol(prelim.data)])
p <- ncol(prelim.data)-1
colnames(prelim.data)[1:p] <- paste("V",1:p,sep="")
colnames(prelim.data)[ncol(prelim.data)]="Y"

#### NO EDITS BELOW THIS LINE (except RESULTS) ####
datamat <- prelim.data

Xmat=datamat[ ,-ncol(datamat)]
Y=datamat[ ,ncol(datamat)]

true.mod=true.model.identification(datamat,thresh = 0.05)

if(runtype=='BMS'){
  res <- est_varselection_bms(datamat = datamat,Xmat = Xmat,true.mod = true.mod,bootn=100)
}else if(runtype=='BMA'){
  res <- est_varselection_study(datamat = datamat,Xmat = Xmat,true.mod = true.mod,bootn=100)
  
}

save(true.mod, file = paste(paste(resultsFile,"truemodel",sep="_"),"rda",sep = "."))

save(res, file = paste(paste(resultsFile,"VSresults",sep="_"),"rda",sep = "."))

write.xlsx(res,file=paste(resultsFile,"xlsx",sep="."))


# data=cbind(model.matrix(Y~.,data=prelim.data),Y)
# data=as.data.frame(data)
# setEPS()
# postscript(paste(paste(resultsFile,"resplot1",sep="_"),"eps",sep = "."))
# par(mfrow=c(1,2))
# plot(true.mod,pch=".",which=c(1,2))
# dev.off()
# 
# setEPS()
# postscript(paste(paste(resultsFile,"resplot2",sep="_"),"eps",sep = "."))
# par(mfrow=c(1,2))
# plot(true.mod,pch=".",which=c(3,5))
# dev.off()
# 
# setEPS()
# postscript(paste(paste(resultsFile,"qqPlot",sep="_"),"eps",sep = "."))
# qqPlot(lm(formula(true.mod),data = data),ylab = "Studentized residuals")
# dev.off()
# 
# 



#load("results/BikeSharingHourly_06252020_truemodel.rda")
# To get tables and plots for true model
#xtable(true.mod)
#stargazer::stargazer(true.mod)

# To load back all the results as list
# library(XLConnect)
# wb <- loadWorkbook("results/BikeSharingHourly_06252020.xlsx")
# res = readWorksheet(wb, sheet = getSheets(wb))

# xtable(res$comp.time,digits=3)
#
# xtable(res$pest.metric*100,digits=3)
#
# xtable(res$uncertainty.coef, digits = 3)
# xtable(res$hypotest.50[ , 1:(ncol(res$hypotest.50)-1)], digits=3)
# xtable(res$hypotest.95[ , 1:(ncol(res$hypotest.95)-1)], digits=3)


