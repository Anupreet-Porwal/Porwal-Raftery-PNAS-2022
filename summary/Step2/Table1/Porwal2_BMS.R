rm(list=ls())
setwd("/mnt/beegfs/homes/porwaa/results/090421/bms/")
load ("VarSel_sum_unnorm_alldata_BMS.rda")
load ("pred_sum_unnorm_alldata_BMS.rda")


R2.res <- pred_sum_unnorm_alldata$R2 
R2.res.norm <-  rowMeans(R2.res)
R2_1.res.norm <- 1-R2.res.norm
R2_1.res.normJZS <- R2_1.res.norm/R2_1.res.norm["JZS-BMA"]


MIS.res <- pred_sum_unnorm_alldata$MeanIntScore.predictions
MIS.res.norm <-  sweep(MIS.res, 2, MIS.res["JZS-BMA" ,],'/')
MISpred.mean <- rowMeans(MIS.res.norm)

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


RMSE.res <- VarSel_sum_unnorm_alldata$RMSE
# Normalise Interval score, width, computation time using JZS performance for each dataset
RMSE.res.norm <-   sweep(RMSE.res, 2, RMSE.res["JZS-BMA", ],'/')
RMSE.mean <- rowMeans(RMSE.res.norm)

IntervalScore.res <- VarSel_sum_unnorm_alldata$ParameterIntScore
IntervalScore.res.norm <-  sweep(IntervalScore.res, 2, IntervalScore.res["JZS-BMA" ,],'/')
MIS.mean <- rowMeans(IntervalScore.res.norm)



AUPRC.res <- VarSel_sum_unnorm_alldata$AUPRC
AUPRC.res.norm <- rowMeans(AUPRC.res)
AUPRC_1.res.norm <- 1-AUPRC.res.norm
AUPRC_1.res.normJZS <- AUPRC_1.res.norm/AUPRC_1.res.norm["JZS-BMA"]



scores <- matrix (rep(NA, 4*5), ncol=5)
rownames (scores) <- rownames(R2.res)
colnames (scores) <- c("RMSE","CI-MIS","AUPRC","1-predR2","predMIS")
scores[,1] <- RMSE.mean # RMSE for point estimates
scores[,2] <- MIS.mean
# MIS for interval estimates
scores[,3] <- AUPRC_1.res.normJZS
scores[,4] <- R2_1.res.normJZS
scores[,5] <- MISpred.mean


agscores <- scores[,1] + scores[,2]+scores[,3]+scores[,4]+scores[,5]
sort (agscores)


# Values are slightly different but should match with BMA values; variability possibly because of sampling 
comptime.res <- VarSel_sum_unnorm_alldata$comptime
comptime.res.norm <-  sweep(comptime.res, 2,comptime.lasso,'/')
comptime.mean <- rowMeans(comptime.res.norm)

phat.res <- pred_sum_unnorm_alldata$phat
phat.res.norm <-  sweep(phat.res, 2, phat.res["JZS-BMA" ,],'/')
phat.mean <- rowMeans(phat.res.norm)


order1 <- c(2,3,1,4)

scores.normJZS.col <- round(
  cbind(agscores[order1]/5,
        scores[order1,],
        phat.mean[order1],
        comptime.mean[order1]),3)
colnames (scores.normJZS.col) <- c("Score","PointEst","IntEst","Inference","Prediction","IntPred","N vars","CPU time")
rownames(scores.normJZS.col) <- c("g=sqrt(n)","Hyper-g","EB-local","JZS-BMA")


datatable (scores.normJZS.col, rownames=TRUE) %>%
  formatStyle(columns = "Score", background = styleInterval(c(0.9999,1.0001, 1.409), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PointEst", background = styleInterval(c(0.9999,1.0001, 1.325), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntEst", background = styleInterval(c(0.9999,1.0001, 1.550), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Inference", background = styleInterval(c(0.9999,1.0001, 1.154), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Prediction", background = styleInterval(c(0.9999,1.0001, 1.238), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntPred", background = styleInterval(c(0.9999,1.0001, 1.174), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "N vars", background = styleInterval(c(0.9999,1.0001, 2.047), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "CPU time", background = styleInterval(c(0.9999,1.0001, 6.678), c("limegreen","lightgrey","yellow","darkorange")))






