setwd("/mnt/beegfs/homes/porwaa/results/")
# 7/12/2021: Producing summary table. See C9236-. Based on porwal.R
load ("VarSel_sum_Norm_alldata.rda")
load ("VarSel_sum_unnorm_alldata.rda")

# RMSE
RMSE.unnorm <- VarSel_sum_unnorm_alldata$RMSE
RMSE.JZS <- RMSE.unnorm["JZS",]
RMSE.normJZS <- RMSE.unnorm %*% diag (1/RMSE.JZS)
colnames (RMSE.normJZS) <- colnames (RMSE.unnorm)
RMSE.normJZSmean <- apply (RMSE.normJZS[-nrow (RMSE.normJZS),], 1, mean)

# MIS for confidence intervals
# In VarSel_sum_Norm_alldata, MIS is already normalized by JZS.
MIS.normJZS <- VarSel_sum_Norm_alldata$ParameterIntScore.Norm
MIS.normJZSmean <- apply (MIS.normJZS[-nrow (MIS.normJZS),], 1, mean)

# AUPRC. 
AUPRC.normTrueModel <- VarSel_sum_unnorm_alldata$AUPRC
AUPRC.normTrueModelmean <- apply (AUPRC.normTrueModel[-nrow (AUPRC.normTrueModel),], 1, mean)

# Now let's normalize it by the JZS value for all the datasets, as for all
#  the other 4 measures. This has the same effect as starting from
#  unnormalized values and then normalizing by JZS.
AUPRC.JZS <- AUPRC.normTrueModel["JZS",]
AUPRC.normJZS <- AUPRC.normTrueModel %*% diag(1/AUPRC.JZS)
colnames (AUPRC.normJZS) <- colnames (AUPRC.normTrueModel)
AUPRC.normJZSmean <- apply (AUPRC.normJZS[-nrow (AUPRC.normJZS),], 1, mean)

# Point predictions: from pred_sum_Norm_alldata. 
load ("pred_sum_Norm_alldata.rda")
predR2 <- pred_sum_Norm_alldata$R2
predR2mean <- apply (predR2[-nrow (predR2),], 1, mean)

# Prediction MIS: MIS is normalized by the JZS value.
predMIS <- pred_sum_Norm_alldata$MeanIntScore.pred.norm
predMISmean <- apply (predMIS[-nrow (predMIS),], 1, mean)

# Form matrix of scores. All scores negatively oriented See C9237.
scores <- matrix (rep(NA, 22*5), ncol=5)
rownames (scores) <- names (RMSE.normJZSmean)
colnames (scores) <- c("RMSE","CI-MIS","AUPRC","1-predR2","predMIS")
scores[,1] <- RMSE.normJZSmean  # RMSE for point estimates
scores[,2] <- c(MIS.normJZSmean[1:12], rep(NA,7), MIS.normJZSmean[13:15])
# MIS for interval estimates
scores[,3] <- c(1-AUPRC.normTrueModelmean[1:17],NA,1-AUPRC.normTrueModelmean[18:21])
# 1 - AUPRC (LASSO-lambda.1se missing for some reason)
scores[,4] <- c(1-predR2mean,NA) # 1 - predR2
scores[,5] <- c(predMISmean[1:12],rep(NA,7),predMISmean[13:14],NA)
# MIS for prediction intervals

# Aggregate scores for methods with intervals
#scores %*% rep(1,5)
agscores <- scores[,1] + scores[,2]+scores[,3]+scores[,4]+scores[,5]
sort (agscores)
plot (agscores[-c(5,9)], rep(0, length(agscores)-2) )

# Aggregate scores for the 3 metrics without intervals
scores[,c(1,3,4)] %*% c(1,1,1)

# Correlation matrix of metrics with only methods that produce intervals
cor (scores[-c(13:19,22),])
pairs (scores[-c(13:19,22),])

cor (scores[-c(18,22),c(1,3,4)])
pairs (scores[-c(18,22),c(1,3,4)])

# 7/8/2021: Normalize the scores matrix so that JZS=1. See C9236-7.
scores.normJZS <- scores
scores.normJZS[,3] <- scores.normJZS[,3] / scores["JZS",3]
scores.normJZS[,4] <- scores.normJZS[,4] / scores["JZS",4]
rownames (scores.normJZS) <- names (RMSE.normJZSmean)
colnames (scores.normJZS) <- c("RMSE","MIS","AUPRC","1-predR2","predMIS")
cor (scores.normJZS[-c(13:19,22),])
pairs (scores.normJZS[-c(13:19,22),])
agscores.normJZS <- scores.normJZS[,1] + scores.normJZS[,2]+scores.normJZS[,3]+scores.normJZS[,4]+scores.normJZS[,5]
sort (agscores.normJZS)

# agscores.normJZS.nonint <- scores.normJZS[,c(1,3,4)] %*% c(1,1,1)
agscores.normJZS.nonint <- scores.normJZS[,1]+scores.normJZS[,3]+scores.normJZS[,4]

# Order the scores matrix in the order on C9238.
cbind (1:22, scores.normJZS)
order1 <- c(8,10,6,21,20,12,3,7,11,17,15,2,19,1,16,13,14,4,9,5,18)
sort (order1) # Check there are no duplicates
round( 
  cbind( 
    agscores.normJZS[order1]/5,
    agscores.normJZS.nonint[order1]/3,
    scores.normJZS[order1,]),
  3)

# Finding thresholds for yellow vs red for the different methods. See C9239.
y <- agscores.normJZS/5
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

y <- agscores.normJZS.nonint/3
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

y <- RMSE.normJZSmean
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

y <- MIS.normJZSmean
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

#y <- 1-AUPRC.normTrueModelmean
y <- scores.normJZS[-22,3]
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

# PredR2
y <- scores.normJZS[-22,4]
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

# PredR2
y <- scores.normJZS[-22,4]
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

# PredR2
y <- scores.normJZS[-22,4]
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

# PredMIS
y <- scores.normJZS[-22,5]
par (mfrow=c(1,2))
plot (1:length(y[!is.na(y)]), sort (y))
plot ((1:(length(y[!is.na(y)])-1)), diff(as.vector(sort (y))))
sort(y)

# Try to produce colored table using DT. 
# https://stackoverflow.com/questions/31323885/how-to-color-specific-cells-in-a-data-frame-table-in-r

#Example:
#library(DT)
#datatable(df, rownames = FALSE) %>%
#  formatStyle(columns = "inputval", 
#              background = styleInterval(c(0.7, 0.8, 0.9)-1e-6, c("white", "lightblue", "magenta", "white"))) %>%
#  formatStyle(columns = "outcome", 
#              background = styleEqual(c(1, 4), c("magenta", "lightblue"))) 

order1 <- c(8,10,6,21,20,12,3,7,11,17,15,2,19,1,16,13,14,4,9,5,18)
scores.normJZS.col <- round(
  cbind(1:21,
        agscores.normJZS[order1]/5,
        agscores.normJZS.nonint[order1]/3,
        scores.normJZS[order1,]),
  3)
colnames (scores.normJZS.col) <- c("Rank","Score","PartScore","RMSE","CI-MIS","AUPRC","1-predR2","pred-MIS")
rownames (scores.normJZS.col) <- 
  c("g=sqrt(n)", "Hyper-g", "EB-local", "ZS-null", "JZS", "Horseshoe", 
    "UIP","EB-global","NLP","LASSO","SCAD","SpikeSlab","Elastic net","BICREG",
    "MCP", "SS Lasso","EMVS","BIC-MCMC","g=1","AIC","LASSO-1se")
scores.normJZS.col

# We use limegreen, yellow and darkorange as the traffic light colors.
#  lightgrey is used to represent the reference JZS method.
library (DT)
datatable (scores.normJZS.col, rownames=TRUE) %>%
  formatStyle(columns = "Score", background = styleInterval(c(0.9999,1.0001, 1.409), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PartScore", background = styleInterval(c(0.9999,1.0001, 1.149), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "RMSE", background = styleInterval(c(0.9999,1.0001, 1.325), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "CI-MIS", background = styleInterval(c(0.9999,1.0001, 1.550), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "AUPRC", background = styleInterval(c(0.9999,1.0001, 1.154), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "1-predR2", background = styleInterval(c(0.9999,1.0001, 1.238), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "pred-MIS", background = styleInterval(c(0.9999,1.0001, 1.174), c("limegreen","lightgrey","yellow","darkorange"))) 

# Need to change order so that Score order is preserved when available.
# Also, AUPRC for Lasso-1se is the same as for Lasso.
# Remove true model
# Remove ZS-null (it's an approx to JZS and gives the same results).
scores.normJZS2 <- cbind (1:22, agscores.normJZS/5, agscores.normJZS.nonint/3,scores.normJZS)
scores.normJZS2 <- scores.normJZS2[-c(nrow(scores.normJZS2)-1,nrow(scores.normJZS2)),] #Remove true model and ZS-null model
colnames (scores.normJZS2) <- c("Rank","Score","PartScore","RMSE","CI-MIS","AUPRC","1-predR2","pred-MIS")
scores.normJZS2["LASSO-lambda.1se","AUPRC"] <- scores.normJZS2["LASSO-lambda.min","AUPRC"]  #Set AUPRC for Lasso-1se is the same as for Lasso.
scores.normJZS2["LASSO-lambda.1se","PartScore"] <- 
  (sum(scores.normJZS2["LASSO-lambda.1se",c(4,6,7)])/3)

# Find revised order.
tmp <- round(cbind (scores.normJZS2[,"Rank"],rank(scores.normJZS2[,"Score"]),
                    rank(scores.normJZS2[,"PartScore"]) ,scores.normJZS2[,"Score"], 
                    scores.normJZS2[,"PartScore"]),3)
colnames (tmp) <-   c("ID","Rank1","Rank2","Score","PartScore")
tmp
tmp[order (tmp[,"Rank1"]),]	# Ordered by score rank (mean of all 5)
tmp[order (tmp[,"Rank2"]),]	# Ordered by non-interval rank (mean of 3)

# Maybe do it by single imputation to imput the full score from the PartScore
y1 <- scores.normJZS2[c(1:12,20),"Score"]
x1 <- scores.normJZS2[c(1:12,20),"PartScore"]
par (mfrow=c(1,1))
plot (x1,y1)
lm1 <- lm (y1~x1)
summary (lm1) #R2=0.79, so correlation = 0.89.
abline (lm1)
par (mfrow=c(2,2))
plot (lm1)

x1missing <- scores.normJZS2[c(13:19),"PartScore"]
y1imp <- -0.1398 + 1.2231 * x1missing 	# Coefficients from lm1
y1imp

# AIC and g=1 are outside the range of the NA methods, so maybe remove them
y2 <- scores.normJZS2[c(1:4,6:8,10:12),"Score"]
x2 <- scores.normJZS2[c(1:4,6:8,10:12),"PartScore"]
par (mfrow=c(1,1))
plot (x2,y2)
lm2 <- lm (y2~x2)
summary (lm2) 
abline (lm2)
par (mfrow=c(2,2))
plot (lm2)
y1imp <- 0.1842 + 0.8795*x1missing	# Coefficients from lm2
y1imp

# New order: see C9238. 
order2 <- c(8,10,6,20,12,3,7,11,17,15,4,1,2,19,16,13,18,14,5,9)
scores.normJZS2.col <- round(cbind(1:nrow(scores.normJZS2), scores.normJZS2[order2,-1]), 3)
scores.normJZS2.col

colnames (scores.normJZS2.col) <- c("Rank","Score","PartScore","PointEst","IntEst","Inference","Prediction","IntPred")
rownames (scores.normJZS2.col) <-
  c("g=sqrt(n)", "Hyper-g", "EB-local", "JZS", "Horseshoe",
    "UIP","EB-global","NLP","LASSO","SCAD","BIC-BAS","BICREG",
    "SpikeSlab","ElasticNet","MCP","SS Lasso","Lasso-1se","EMVS","AIC","g=1")
scores.normJZS2.col

# We use limegreen, yellow and darkorange as the traffic light colors.
#  lightgrey is used to represent the reference JZS method.
library (DT)
datatable (scores.normJZS2.col, rownames=TRUE) %>%
  formatStyle(columns = "Score", background = styleInterval(c(0.9999,1.0001, 1.409), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PartScore", background = styleInterval(c(0.9999,1.0001, 1.148), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PointEst", background = styleInterval(c(0.9999,1.0001, 1.325), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntEst", background = styleInterval(c(0.9999,1.0001, 1.550), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Inference", background = styleInterval(c(0.9999,1.0001, 1.154), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Prediction", background = styleInterval(c(0.9999,1.0001, 1.238), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntPred", background = styleInterval(c(0.9999,1.0001, 1.174), c("limegreen","lightgrey","yellow","darkorange")))

# 7/16/2021: Add computer time to the results
VarSel_sum_Norm_alldata$comptime
comptime.unnorm <- VarSel_sum_unnorm_alldata$comptime
comptime.LASSO <- comptime.unnorm["LASSO-lambda.min",]
comptime.normLASSO <- comptime.unnorm %*% diag(1/comptime.LASSO)
comptime.normLASSOmean <- apply(comptime.normLASSO, 1, mean)[1:20]
comptime.normLASSOmean <- comptime.normLASSOmean[order2]
comptime.normLASSOmean

# Determine breakpoints for colors
par (mfrow=c(1,2))
plot ((1:length(comptime.normLASSOmean)), sort (comptime.normLASSOmean))
plot ((1:(length(comptime.normLASSOmean)-1)), diff(as.vector(sort (comptime.normLASSOmean))))
sort (comptime.normLASSOmean)

# Add phat to the results
load ("pred_sum_Norm_alldata.rda")
phat <- pred_sum_Norm_alldata$phat
phat.JZS <- phat["JZS",]
phat.normJZS <- phat %*% diag (1/phat.JZS)
phat.normJZS
phat.normJZSmean <- apply (phat.normJZS,1,mean)
phat.normJZSmean <- phat.normJZSmean[1:20]
phat.normJZSmean <- phat.normJZSmean[order2]
phat.normJZSmean

# Determine breakpoints for colors
par (mfrow=c(1,2))
plot ((1:length(phat.normJZSmean)), sort (phat.normJZSmean))
plot ((1:(length(phat.normJZSmean)-1)), diff(as.vector(sort (phat.normJZSmean))))
sort (phat.normJZSmean)

order2 <- c(8,10,6,20,12,3,7,11,17,15,4,1,2,19,16,13,18,14,5,9)
scores.normJZS3.col <- round(cbind(1:nrow(scores.normJZS2), scores.normJZS2[order2,-1],phat.normJZSmean,comptime.normLASSOmean), 3)
scores.normJZS3.col

colnames (scores.normJZS3.col) <- c("Rank","Score","PartScore","PointEst","IntEst","Inference","Prediction","IntPred","N vars","CPU time")
rownames (scores.normJZS3.col) <-
  c("g=sqrt(n)", "Hyper-g", "EB-local", "JZS", "Horseshoe",
    "UIP","EB-global","NLP","LASSO","SCAD","BIC-BAS","BICREG",
    "SpikeSlab","ElasticNet","MCP","SS Lasso","Lasso-1se","EMVS","AIC","g=1")
scores.normJZS3.col

library (DT)
result <-  datatable (scores.normJZS3.col, rownames=TRUE) %>%
  formatStyle(columns = "Score", background = styleInterval(c(0.9999,1.0001, 1.409), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PartScore", background = styleInterval(c(0.9999,1.0001, 1.148), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "PointEst", background = styleInterval(c(0.9999,1.0001, 1.325), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntEst", background = styleInterval(c(0.9999,1.0001, 1.550), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Inference", background = styleInterval(c(0.9999,1.0001, 1.154), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "Prediction", background = styleInterval(c(0.9999,1.0001, 1.238), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "IntPred", background = styleInterval(c(0.9999,1.0001, 1.174), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "N vars", background = styleInterval(c(0.9999,1.0001, 2.047), c("limegreen","lightgrey","yellow","darkorange"))) %>%
  formatStyle(columns = "CPU time", background = styleInterval(c(0.9999,1.0001, 6.678), c("limegreen","lightgrey","yellow","darkorange")))

htmlwidgets::saveWidget(result, "result.html")