rm(list=ls())
setwd("/mnt/beegfs/homes/porwaa/results/090421/bms/")


dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily", 
               "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston", 
               "Nutrimouse","Multidrug","NIR","Liver")

resultsFile = paste(dataname,"_","prediction","_","BMS",sep = "")

filename <- paste(resultsFile,"rda",sep = ".")


methods_Point=c("EB-local","g-sqrtn","Hyper g","JZS-BMA")# ,"ScanBMA")
methods_Interval=c("EB-local","g-sqrtn","Hyper g","JZS-BMA")

R2.res <- matrix(NA, length(methods_Point),length(dataname))
rownames(R2.res) <- methods_Point
colnames(R2.res) <- dataname


phat.res <- matrix(NA, length(methods_Point),length(dataname))
rownames(phat.res) <- methods_Point
colnames(phat.res) <- dataname

Coverage.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(Coverage.res) <- methods_Interval
colnames(Coverage.res) <- dataname

MIS.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(MIS.res) <- methods_Interval
colnames(MIS.res) <- dataname

Width.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(Width.res) <- methods_Interval
colnames(Width.res) <- dataname


CRPS.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(CRPS.res) <- methods_Interval
colnames(CRPS.res) <- dataname

for (i in 1:length(filename)){
  load(file=filename[i])
  R2.res[ ,i] <- res1$point.pred.metric[ ,4]
  phat.res[ ,i] <- res1$point.pred.metric[ ,5]
  Coverage.res[ ,i] <- res1$uncertainty.pred.metric[ ,1]
  Width.res[ ,i] <- res1$uncertainty.pred.metric[ ,2]
  MIS.res[ ,i] <- res1$uncertainty.pred.metric[ ,3]
  CRPS.res[ ,i] <- res1$uncertainty.pred.metric[ ,4]
}

R2.res.norm <-  rowMeans(R2.res)
R2_1.res.norm <- 1-R2.res.norm
R2_1.res.normJZS <- R2_1.res.norm/R2_1.res.norm["JZS-BMA"]




phat.res.norm <-  sweep(phat.res, 2, phat.res["JZS-BMA" ,],'/')

#Coverage.error.res <- abs(Coverage.res -0.95)
#Width.res.norm <-  sweep(Width.res, 2, Width.res[4,],'/')
MIS.res.norm <-  sweep(MIS.res, 2, MIS.res["JZS-BMA" ,],'/')
#CRPS.res.norm <-  sweep(CRPS.res, 2, CRPS.res[4 ,],'/')


pred_sum_unnorm_alldata <-  list("R2"=R2.res, "phat"=phat.res, "PredCoverage"=Coverage.res,
                                 "MeanIntScore.predictions"=MIS.res, "CRPS.predictions"=CRPS.res,"Width.predictions"=Width.res)


round(R2_1.res.normJZS[c(2,3,1,4)],digits=3)
round(rowMeans(MIS.res.norm)[c(2,3,1,4)],digits=3)
round(rowMeans(phat.res.norm)[c(2,3,1,4)],digits=3)



# pred_sum_Norm_alldata <-  list("R2"=R2.res.norm, 
#                                "phat"=phat.res.norm, 
#                                #"PredCoverageerror"=Coverage.error.res, 
#                                "MeanIntScore.pred.norm"=MIS.res.norm#, 
#                                #"CRPS.pred.norm"=CRPS.res.norm,
#                                #"Width.pred.norm"=Width.res.norm
#                                )
# 
# 
# 
# 
save(pred_sum_unnorm_alldata, file = "pred_sum_unnorm_alldata_BMS.rda")
# 
# save(pred_sum_Norm_alldata, file = "pred_sum_Norm_alldata_BMS.rda")
# 
# 
