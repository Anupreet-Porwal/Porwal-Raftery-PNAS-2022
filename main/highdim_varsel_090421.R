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
data_with_int_quad <- function(original.data){
  # Standardizing all the x's except the factor variables
  for (i in 1:(ncol(original.data)-1)){
    if(is.factor(original.data[ ,i])==0){
      original.data[ ,i]=(original.data[ ,i]-mean(original.data[ ,i]))/sd(original.data[ ,i])
    }
  }

  # centering the Y. Assuming Y is always the last column of original.data
  y.val.0 <- mean(original.data[ , ncol(original.data)])
  original.data[ , ncol(original.data)]=original.data[ , ncol(original.data)]-mean(original.data[ , ncol(original.data)])
  p <- ncol(original.data)-1
  colnames(original.data)[ncol(original.data)]="Y"
  colnames(original.data)[1:p] <- paste("V",1:p,sep="")


  X.original <- original.data[ ,-ncol(original.data)]
  Y.original <- original.data[ , ncol(original.data)]

  colnum <- ncol(X.original)
  coltype <- sapply(X.original, class)


  # Include all the two way interaction terms; remove the intercept term in next line
  X.original <- model.matrix(Y~.*.,data=original.data)
  X.original <- X.original[ ,-1]
  # for (j in 1:(colnum-1)){
  #   for (k in (j+1):colnum){
  #     X.original <- cbind.data.frame(X.original, X.original[ ,j]*X.original[ ,k])
  #     colnames(X.original)[ncol(X.original)] <-
  #       paste(colnames(X.original)[j],colnames(X.original)[k],sep=".")
  #   }
  #
  # }
  #
  for (j in 1:colnum){
    if(coltype[j]=="numeric"){
      X.original <- cbind.data.frame(X.original, X.original[ ,j]^2)
      colnames(X.original)[ncol(X.original)] <-
        paste(colnames(X.original)[j],2,sep=".")
    }else if(coltype[j]=="factor"){
      X.original[ ,j] <- as.factor(X.original[ ,j])
    }
  }


  prelim.data <- cbind.data.frame(X.original,Y.original)
  colnames(prelim.data)[ncol(prelim.data)] <- "Y"

  return(prelim.data)

}

#setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")
#source("src/est_varselection_function_SIS_011321.R")

args = commandArgs(TRUE)

data.num <- as.numeric(args[1])
runtype <- as.character(args[3])


if(runtype=='BMS'){
  source("src/bms_est_varselection_function_SIS_083121.R")
}else if(runtype=='BMA'){
  source("src/est_varselection_function_SIS_011321.R")
}


datanames <- c("Diabetes", "Superconductivity", "Ozone", "Boston",
               "Nutrimouse","Multidrug","NIR","Liver")


dataname <- datanames[data.num]


resultsFile = paste(as.character(args[2]), paste0(dataname,
                                                  "_", runtype),sep = "")



if(dataname=="Diabetes"){
  data("diabetesI")
  prelim.data <- diabetesI

  prelim.data <- cbind.data.frame(prelim.data[ ,-1],prelim.data[ ,1])

}else if(dataname=="Superconductivity"){
  prelim.data <- read.table("data/Superconductivity/train.csv",header = TRUE,sep = ",")
  # transforming Y variable with cube root transformation
  prelim.data$critical_temp <- prelim.data$critical_temp^(1/3)

}else if(dataname=="Ozone"){
  # Data version used in Liang et Al (2008), Miller 2001
  data("ozone",package = "gss")
  ozone <- ozone[ ,-ncol(ozone)]
  original.data <- cbind.data.frame(ozone[ ,-1], ozone[ ,1])

  colnames(original.data)[ncol(original.data)] <- "Ozone"
  original.data$Ozone <- log(original.data$Ozone)

  prelim.data <- data_with_int_quad(original.data)



}else if(dataname=="Boston"){
  data("BostonHousing")
  original.data <- BostonHousing
  prelim.data <- data_with_int_quad(original.data)

}else if(dataname=="Nutrimouse"){
  library(mixOmics)

  data("nutrimouse")

  prelim.data <- cbind(nutrimouse$gene,nutrimouse$lipid$C16.0)
  prelim.data <- as.data.frame(prelim.data)
  colnames(prelim.data)[ncol(prelim.data)] <- colnames(nutrimouse$lipid)[2]

}else if(dataname=="Multidrug"){
  library(mixOmics)
  data("multidrug")

  x.dirty <- multidrug$compound
  x.clean <- x.dirty[ , colSums(is.na(x.dirty))==0]

  y <- multidrug$ABC.trans[ ,3]

  prelim.data <- cbind(x.clean,y)
  prelim.data <- as.data.frame(prelim.data)

}else if(dataname=="NIR"){
  library(chemometrics)

  data(NIR)
  # response is square root of glucose
  prelim.data <- cbind(NIR$xNIR,sqrt(NIR$yGlcEtOH$Glucose))
  prelim.data <- as.data.frame(prelim.data)

}else if(dataname=="Liver"){
  library(mixOmics)
  data("liver.toxicity")

  prelim.data <- cbind(liver.toxicity$gene, liver.toxicity$clinic$Cholesterol.mg.dL.)
  prelim.data <- as.data.frame(prelim.data)

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



