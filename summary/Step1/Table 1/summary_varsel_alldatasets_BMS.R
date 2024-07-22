rm(list=ls())

setwd("/mnt/beegfs/homes/porwaa/results/090421/bms/")

dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily", 
              "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston", 
              "Nutrimouse","Multidrug","NIR","Liver")


resultsFile = paste(dataname,'_','BMS',sep = "")
filename <- paste(paste(resultsFile,"VSresults",sep="_"),"rda",sep = ".")


methods=c("EB-local","g-sqrtn","Hyper g","JZS-BMA")


AUPRC.res <- matrix(NA, length(methods),length(dataname))
rownames(AUPRC.res) <- methods
colnames(AUPRC.res) <- dataname

AUROC.res <- matrix(NA, length(methods),length(dataname))
rownames(AUROC.res) <- methods
colnames(AUROC.res) <- dataname

F1.res <- matrix(NA, length(methods),length(dataname))
rownames(F1.res) <- methods
colnames(F1.res) <- dataname

TER.res <- matrix(NA, length(methods),length(dataname))
rownames(TER.res) <- methods
colnames(TER.res) <- dataname

RMSE.res <- matrix(NA, length(methods),length(dataname))
rownames(RMSE.res) <- methods
colnames(RMSE.res) <- dataname

oracle.rmse <- matrix(NA, 1,length(dataname))
colnames(oracle.rmse) <- dataname


comptime.res <- matrix(NA, length(methods),length(dataname))
rownames(comptime.res) <- methods
colnames(comptime.res) <- dataname


Coverage.res <- matrix(NA, length(methods),length(dataname))
rownames(Coverage.res) <- methods
colnames(Coverage.res) <- dataname

zeroCoverage.res <- matrix(NA, length(methods),length(dataname))
rownames(zeroCoverage.res) <- methods
colnames(zeroCoverage.res) <- dataname

nzeroCoverage.res <- matrix(NA, length(methods),length(dataname))
rownames(nzeroCoverage.res) <- methods
colnames(nzeroCoverage.res) <- dataname


Width.res <- matrix(NA, length(methods),length(dataname))
rownames(Width.res) <- methods
colnames(Width.res) <- dataname


IntervalScore.res <- matrix(NA, length(methods),length(dataname))
rownames(IntervalScore.res) <- methods
colnames(IntervalScore.res) <- dataname



for (i in 1:length(dataname)){
  load(file = filename[i])
  AUPRC.res[ ,i]=colMeans(res$AUPRC.mat)
  AUROC.res[ ,i]=colMeans(res$AUROC.mat)
  F1.res[ ,i]= res$hypotest[ -nrow(res$hypotest) ,4]
  TER.res[ ,i]= res$hypotest[ -nrow(res$hypotest) , 7]
  # load(file=truemod.filename[i])
  # oracle.rmse[ ,i] <- sqrt(mean(coef(summary(true.mod))[-1 ,2]^2))
  RMSE.res[ ,i]=res$pest.metric[ ,3]
  comptime.res[ ,i]=res$comp.time
  Coverage.res[ ,i]=res$uncertainty.coef[1, ]
  Width.res[ ,i] = res$uncertainty.coef[2, ]
  IntervalScore.res[ ,i]=res$uncertainty.coef[3, ]
  zeroCoverage.res[ ,i]=res$zero_coverage
  nzeroCoverage.res[ ,i]=res$nzero_coverage
}

# taken from Varsel_sum_unnorm_alldata$comptime
comptime.lasso <- c(0.11064,
0.27689,
0.27899,
0.14369,
0.46994,
0.12323,
0.33156,
1.67711,
0.16260,
0.50055,
0.12053,
0.20861,
0.40362,
0.39693)


# Normalise Interval score, width, computation time using JZS performance for each dataset
RMSE.res.norm <-   sweep(RMSE.res, 2, RMSE.res["JZS-BMA", ],'/')
IntervalScore.res.norm <-  sweep(IntervalScore.res, 2, IntervalScore.res["JZS-BMA" ,],'/')
comptime.res.norm <-  sweep(comptime.res, 2,comptime.lasso,'/')
#Width.res.norm <-  sweep(Width.res, 2, Width.res[4,],'/')
#coverage.error.res <- abs(Coverage.res-0.95)

AUPRC.res.norm <- rowMeans(AUPRC.res)
AUPRC_1.res.norm <- 1-AUPRC.res.norm
AUPRC_1.res.normJZS <- AUPRC_1.res.norm/AUPRC_1.res.norm["JZS-BMA"]


round(rowMeans(RMSE.res.norm)[c(2,3,1,4)],digits=3)
round(rowMeans(IntervalScore.res.norm)[c(2,3,1,4)],digits=3)
round(AUPRC_1.res.normJZS[c(2,3,1,4)],digits=3)
# Comptime for BMS is same as BMA since we use the same function
#round(rowMeans(comptime.res.norm)[c(2,3,1,4)],digits=3)



VarSel_sum_unnorm_alldata <- list("RMSE"=RMSE.res,
                                  "oracleRMSE"=oracle.rmse , "ParameterCoverage"=Coverage.res,
                                  "ParameterIntScore"=IntervalScore.res, "ParameterWidth"=Width.res,
                                  "AUPRC"=AUPRC.res, "AUROC"=AUROC.res
                                  ,"comptime"=comptime.res, "F1score"=F1.res, "TER"=TER.res)
# 
# VarSel_sum_Norm_alldata <- list("RMSE.Norm"=RMSE.res.norm,
#                                 #"oracleRMSE"=oracle.rmse , 
#                                 #"ParameterCoverageError"=coverage.error.res, 
#                                 "ParameterIntScore.Norm"=IntervalScore.res.norm, 
#                                 #"ParameterWidth.Norm"=Width.res.norm,
#                                 "AUPRC.Norm"=AUPRC.res.norm, 
#                                 #"AUROC.Norm"=auroc.res.normalised,
#                                 "comptime.Norm"=comptime.res.norm)



save(VarSel_sum_unnorm_alldata, file = "VarSel_sum_unnorm_alldata_BMS.rda")
# save(VarSel_sum_Norm_alldata, file = "VarSel_sum_Norm_alldata.rda")
