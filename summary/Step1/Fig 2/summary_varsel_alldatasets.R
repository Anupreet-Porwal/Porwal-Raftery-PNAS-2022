setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")

dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily", 
               "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston", 
               "Nutrimouse","Multidrug","NIR","Liver")


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



resultsFile = paste("results/113021/bma/", paste0(dataname, 
                                       "_", "BMA"),sep = "")
filename <- paste(paste(resultsFile,"VSresults",sep="_"),"rda",sep = ".")

truemod.file =paste("results/113021/bma/", paste0(dataname, 
                                              "_", "BMA"),sep = "")
truemod.filename <- paste(paste(truemod.file,"truemodel",sep="_"),"rda",sep = ".")

methods_varsel_prc= c("BMA-bicreg","Spikeslab", "UIP", "Benchmark","BIC","AIC","EB-local","EB-global",
                      "g-sqrt n", "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS","SCAD",
                      "MCP","LASSO","Elastic net", "JZS","ZS-null", "true.mod" , "full model")

methods_varsel=c("BMA-bicreg","Spikeslab", "UIP", "Benchmark","BIC","AIC","EB-local","EB-global","g-sqrt n",
                 "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS","SCAD","MCP",
                 "LASSO-lambda.min","LASSO-lambda.1se","Elastic net", "JZS","ZS-null", 
                 "true.mod" , "full model")


methods_Point=c("BMA-bicreg","Spikeslab", "UIP", "Benchmark", "BIC","AIC","EB-local","EB-global","g-sqrt n",
                "g-1","Hyper g","NLP","Horseshoe","SS Lasso","EMVS","SCAD","MCP","LASSO-lambda.min",
                "LASSO-lambda.1se","Elastic net", "JZS","ZS-null","true.mod", "full model")

methods_Interval=c("BMA -bicreg","Spikeslab","UIP", "Benchmark", "BIC", "AIC","EB-local","EB-global",
                   "g-sqrt n" , "g-1","Hyper g", "NLP", "Horseshoe",  "JZS","ZS-null" , 
                   "true.mod", "full model")#,"ScanBMA")


AUPRC.res <- matrix(NA, length(methods_varsel_prc),length(dataname))
rownames(AUPRC.res) <- methods_varsel_prc
colnames(AUPRC.res) <- dataname

AUROC.res <- matrix(NA, length(methods_varsel_prc),length(dataname))
rownames(AUROC.res) <- methods_varsel_prc
colnames(AUROC.res) <- dataname

F1.res <- matrix(NA, length(methods_varsel),length(dataname))
rownames(F1.res) <- methods_varsel
colnames(F1.res) <- dataname

TER.res <- matrix(NA, length(methods_varsel),length(dataname))
rownames(TER.res) <- methods_varsel
colnames(TER.res) <- dataname


RMSE.res <- matrix(NA, length(methods_Point),length(dataname))
rownames(RMSE.res) <- methods_Point
colnames(RMSE.res) <- dataname

oracle.rmse <- matrix(NA, 1,length(dataname))
colnames(oracle.rmse) <- dataname


comptime.res <- matrix(NA, length(methods_Point),length(dataname))
rownames(comptime.res) <- methods_Point
colnames(comptime.res) <- dataname


Coverage.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(Coverage.res) <- methods_Interval
colnames(Coverage.res) <- dataname

zeroCoverage.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(zeroCoverage.res) <- methods_Interval
colnames(zeroCoverage.res) <- dataname

nzeroCoverage.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(nzeroCoverage.res) <- methods_Interval
colnames(nzeroCoverage.res) <- dataname


Width.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(Width.res) <- methods_Interval
colnames(Width.res) <- dataname


IntervalScore.res <- matrix(NA, length(methods_Interval),length(dataname))
rownames(IntervalScore.res) <- methods_Interval
colnames(IntervalScore.res) <- dataname



for (i in 1:length(dataname)){
  load(file = filename[i])
  AUPRC.res[ ,i]=colMeans(res$AUPRC.mat)
  AUROC.res[ ,i]=colMeans(res$AUROC.mat)
  F1.res[ ,i]= res$hypotest.50[ -nrow(res$hypotest.50) ,4]
  TER.res[ ,i]= res$hypotest.50[ -nrow(res$hypotest.50) , 7]
  load(file=truemod.filename[i])
  oracle.rmse[ ,i] <- sqrt(mean(coef(summary(true.mod))[-1 ,2]^2))
  RMSE.res[ ,i]=res$pest.metric[ ,3]#/oracle.rmse[ ,i]
  comptime.res[ ,i]=res$comp.time
  Coverage.res[ ,i]=res$uncertainty.coef[1, ]
  Width.res[ ,i] = res$uncertainty.coef[2, ]
  IntervalScore.res[ ,i]=res$uncertainty.coef[3, ]
  zeroCoverage.res[ ,i]=res$zero_coverage
  nzeroCoverage.res[ ,i]=res$nzero_coverage
}

# AUPRC baseline is given by P/(P+N) with in our case becomes p0/p
# Soruce: https://stats.stackexchange.com/questions/251175/what-is-baseline-in-precision-recall-curve
auprc.baseline <- dat.dim[ ,3]/dat.dim[ ,2]

auprc.res.normalised <- sweep(AUPRC.res, 2, auprc.baseline, `/`)

# This is fixed at 0.5 irrespective of the dataset
auroc.baseline <- 0.5

auroc.res.normalised <- AUROC.res/auroc.baseline

# # F1 baseline is given as follows if s=p0/p then s/(s+0.5) is the F1 score for random classifier
# # classifying one with probability 0.5 
# # Source: https://stats.stackexchange.com/questions/390200/what-is-the-baseline-of-the-f1-score-for-a-binary-classifier
# s <- dat.dim[ ,3]/dat.dim[ ,2]
# F1.baseline <- s/(s+0.5)
# F1.res.normalised <- sweep(F1.res, 2, F1.baseline, `/`)

# Normalise Interval score, width, computation time using JZS performance for each dataset
IntervalScore.res.norm <-  sweep(IntervalScore.res, 2, IntervalScore.res[14 ,],'/')
Width.res.norm <-  sweep(Width.res, 2, Width.res[ 14,],'/')
RMSE.res.norm <-   sweep(RMSE.res, 2, oracle.rmse,'/')
comptime.res.norm <-  sweep(comptime.res, 2, comptime.res[ 21,],'/')
coverage.error.res <- abs(Coverage.res-0.95)



VarSel_sum_unnorm_alldata <- list("dataset.true.mod.summary" =dat.dim,"RMSE"=RMSE.res,
                         "oracleRMSE"=oracle.rmse , "ParameterCoverage"=Coverage.res, 
                         "ParameterIntScore"=IntervalScore.res, "ParameterWidth"=Width.res,
                         "AUPRC"=AUPRC.res, "AUROC"=AUROC.res
                         ,"comptime"=comptime.res, "F1score"=F1.res, "TER"=TER.res)

VarSel_sum_Norm_alldata <- list("RMSE.Norm"=RMSE.res.norm,
            "oracleRMSE"=oracle.rmse , "ParameterCoverageError"=coverage.error.res, 
            "ParameterIntScore.Norm"=IntervalScore.res.norm, "ParameterWidth.Norm"=Width.res.norm,
            "AUPRC.Norm"=auprc.res.normalised, "AUROC.Norm"=auroc.res.normalised,
            "comptime.Norm"=comptime.res.norm)


save(VarSel_sum_unnorm_alldata, file = "Summary-Rcode/VarSel_sum_unnorm_alldata.rda")
save(VarSel_sum_Norm_alldata, file = "Summary-Rcode/VarSel_sum_Norm_alldata.rda")
