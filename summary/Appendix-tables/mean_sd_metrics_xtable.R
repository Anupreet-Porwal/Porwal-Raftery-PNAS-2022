library(xtable)


setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")

dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily", 
              "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston", 
              "Nutrimouse","Multidrug","NIR","Liver")
#dates_vars <- rep("01142021", length(dataname))


dat.dim <- matrix(NA, nrow=length(dataname),ncol = 4)
rownames(dat.dim) <- dataname
colnames(dat.dim) <- c("n","p","p0","R2")

dat.dim[1, ] <- c(777,14,7,0.8659)
dat.dim[2, ] <- c(7590,21,17,0.7809)
dat.dim[3, ] <- c(7590,21,20,0.8375)
dat.dim[4, ] <- c(731,28,16,0.8439)
dat.dim[5, ] <- c(17379,32,21,0.7238)
dat.dim[6, ] <- c(1373,22,16,0.9306)
dat.dim[7, ] <- c(442,64,7,0.5171)
dat.dim[8, ] <- c(21263,81,27,0.7216)
dat.dim[9, ] <- c(330,44,5,0.7489)
dat.dim[10, ] <- c(506,103,23,0.8621)
dat.dim[11, ] <- c(40,120,9,0.9447)
dat.dim[12, ] <- c(60,853,12,0.8490)
dat.dim[13, ] <- c(166,225,7,0.8648)
dat.dim[14, ] <- c(64,3116,8,0.844)


mean_sd_function <- function(filename,dataname, methods.names, metric,mag.factor,dig,neg=FALSE){
  res.mat <- matrix(NA, length(methods.names),4*length(dataname))
  rownames(res.mat) <- methods.names
  
  for (j in 1:length(dataname)){
    ind <- (j-1)*4
    load(file = filename[j])
    if(neg==TRUE){
      res.metric <- 1-res[[metric]]
    }else{
      res.metric <- res[[metric]]
    }
    res.mat[ ,ind+1] <- round(apply(res.metric , 2, mean)*mag.factor,digits=dig)
    res.mat[ ,ind+2] <- "bop" 
    res.mat[ ,ind+3] <- round(apply(res.metric, 2, sd)*mag.factor,digits=dig)
    res.mat[ ,ind+4] <-  "bcl"
    
}
  return(res.mat)
  
}
resultsFile = paste("results/113021/bma/", paste0(dataname, 
                                              "_","BMA"),sep = "")
filename <- paste(paste(resultsFile,"VSresults",sep="_"),"rda",sep = ".")

truemod.file =paste("results/113021/bma/", paste0(dataname, 
                                              "_", "BMA"),sep = "")
truemod.filename <- paste(paste(truemod.file,"truemodel",sep="_"),"rda",sep = ".")

methods_varsel_prc= c("BMA-bicreg","Spikeslab", "UIP", "Benchmark", "BIC","AIC","EB-local","EB-global",
                      "g-sqrt n", "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS","SCAD",
                      "MCP","LASSO","Elastic net", "JZS","ZS-null", "true.mod" , "full model")

methods_varsel=c("BMA-bicreg","Spikeslab", "UIP", "Benchmark", "BIC","AIC","EB-local","EB-global","g-sqrt n",
                 "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS","SCAD","MCP",
                 "LASSO-lambda.min","LASSO-lambda.1se","Elastic net", "JZS","ZS-null", 
                 "true.mod" , "full model")


methods_Point=c("BMA-bicreg","Spikeslab", "UIP", "Benchmark", "BIC","AIC","EB-local","EB-global","g-sqrt n",
                "g-1","Hyper g","NLP","Horseshoe","SS Lasso","EMVS","SCAD","MCP","LASSO-lambda.min",
                "LASSO-lambda.1se","Elastic net", "JZS","ZS-null","true.mod", "full model")

methods_Interval=c("BMA -bicreg","Spikeslab","UIP", "Benchmark", "BIC", "AIC","EB-local","EB-global",
                   "g-sqrt n" , "g-1","Hyper g", "NLP", "Horseshoe",  "JZS","ZS-null" , 
                   "true.mod", "full model")#,"ScanBMA")

# Calculating 1-AUPRC and 1-AUROC
AUPRC.res <- mean_sd_function(filename,dataname,methods.names=methods_varsel_prc,metric = "AUPRC.mat",mag.factor = 1,dig = 3,neg = TRUE)
AUROC.res <- mean_sd_function(filename,dataname,methods.names=methods_varsel_prc,metric = "AUROC.mat",mag.factor = 1,dig = 3,neg = TRUE)

F1.res <- mean_sd_function(filename,dataname,methods.names=methods_varsel,metric = "F1.mat5",mag.factor = 1,dig = 4)
TER.res <- mean_sd_function(filename,dataname,methods.names=methods_varsel,metric = "total.err.mat5",mag.factor = 1,dig = 4)


RMSE.res <- mean_sd_function(filename,dataname,methods.names=methods_Point,metric = "rmse_mat_coef",mag.factor = 1000,dig = 3)
comptime.res <- mean_sd_function(filename,dataname,methods.names=methods_Point,metric = "comp_time_fit",mag.factor = 1,dig = 3)


coverage.res <- mean_sd_function(filename,dataname,methods.names=methods_Interval,metric = "coverage_mat_coef",mag.factor = 1,dig = 4)
IntervalScore.res <- mean_sd_function(filename,dataname,methods.names=methods_Interval,metric = "IntervalScore_mat_coef",mag.factor = 1,dig = 3)

order.all <- c(9,11,7,21,13,3,8,12,4,18,16,5,1,2,20,17,14,19,15,6,10)

order.int <- c(9,11,7,14,13,3,8,12,4,5,1,2,6,10)

order.prc <- c(9,11,7,20,13,3,8,12,4,18,16,5,1,2,19,17,14,15,6,10)

xtable(RMSE.res[order.all, ])
xtable(coverage.res)
xtable(IntervalScore.res[order.int, ])
xtable(AUPRC.res[order.prc, ])
xtable(F1.res[ , 1:40])
xtable(TER.res[ ,1:40])
xtable(comptime.res[order.all, ])


rm(list=ls())

#____________________________________
dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily","BS-hourly","SML","Diabetes","Superconductivity",
              "Ozone", "Boston", "Nutrimouse","Multidrug","NIR","Liver") 
# dates_pred <- rep("01142021", length(dataname))
# dates_pred[5] <- "01202021"
# dates_pred[8] <- "01232021"

resultsFile = paste("results/113021/bma/", paste0(dataname,"_","prediction","_", 
                                              "BMA"),sep = "")

filename <- paste(resultsFile,"rda",sep = ".")


methods_Point=c("BMA-bicreg","Spikeslab","UIP","Benchmark", "BIC ","AIC","EB-local","EB-global",
                "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","SS Lasso", "EMVS", "SCAD",
                "MCP","LASSO-lambda.min","LASSO-lambda.1se","Elastic Net","JZS","ZS-null","full model")# ,"ScanBMA")
methods_Interval=c("BMA-bicreg","Spikeslab","UIP","Benchmark", "BIC ","AIC","EB-local","EB-global",
                   "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","JZS","ZS-null","full model")# ,"ScanBMA")


mean_sd_function_pred <- function(filename,dataname, methods.names, metric,mag.factor,dig,neg=FALSE){
  res.mat <- matrix(NA, length(methods.names),4*length(dataname))
  rownames(res.mat) <- methods.names
  
  for (j in 1:length(dataname)){
    ind <- (j-1)*4
    load(file = filename[j])
    if(neg==TRUE){
      res.metric <- 1-res1[[metric]]
    }else{
      res.metric <- res1[[metric]]
    }
    
    res.mat[ ,ind+1] <- round(apply(res.metric, 2, mean)*mag.factor,digits=dig)
    res.mat[ ,ind+2] <- "bop" 
    res.mat[ ,ind+3] <- round(apply(res.metric, 2, sd)*mag.factor,digits=dig)
    res.mat[ ,ind+4] <-  "bcl"
    
  }
  return(res.mat)
  
}

# Calculate 1-R2
R2.res <- mean_sd_function_pred(filename,dataname,methods_Point,metric = "r2test_mat",mag.factor = 1,dig=3,neg=TRUE)
phat.res <- mean_sd_function_pred(filename,dataname,methods_Point,metric = "phat.mat",mag.factor = 1,dig=2)
Coverage.res <- mean_sd_function_pred(filename,dataname,methods_Interval,metric = "coverage_mat",mag.factor = 1,dig=3)
MIS.res <- mean_sd_function_pred(filename,dataname,methods_Interval,metric = "IntervalScore",mag.factor = 1,dig=3)
#CRPS.res <- mean_sd_function_pred(filename,dataname,methods_Interval,metric = "CRPS_mat",mag.factor = 1,dig=3)

order.all <- c(9,11,7,21,13,3,8,12,4,18,16,5,1,2,20,17,14,19,15,6,10)

order.int <- c(9,11,7,14,13,3,8,12,4,5,1,2,6,10)

order.prc <- c(9,11,7,20,13,3,8,12,4,18,16,5,1,2,19,17,14,15,6,10)

xtable(R2.res[order.all, ])
xtable(phat.res[order.all, ])
xtable(MIS.res[order.int, ])
