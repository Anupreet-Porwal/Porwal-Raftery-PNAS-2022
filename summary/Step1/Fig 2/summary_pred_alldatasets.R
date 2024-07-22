setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")


datanames <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily", 
               "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston", 
               "Nutrimouse","Multidrug","NIR","Liver")

dat.dim <- matrix(NA, nrow=length(datanames),ncol = 4)
rownames(dat.dim) <- datanames
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




dataname <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily","BS-hourly","SML","Diabetes","Superconductivity",
               "Ozone", "Boston", "Nutrimouse","Multidrug","NIR","Liver") 
# dates_pred <- rep("01142021", length(dataname))
# dates_pred[5] <- "01202021"
# dates_pred[8] <- "01232021"

resultsFile = paste("results/113021/bma/", paste0(dataname,"_","prediction","_", 
                              "BMA"),sep = "")

filename <- paste(resultsFile,"rda",sep = ".")


methods_Point=c("BMA-bicreg","Spikeslab","UIP", "Benchmark", "BIC ","AIC","EB-local","EB-global",
                "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","SS Lasso", "EMVS", "SCAD",
                "MCP","LASSO-lambda.min","LASSO-lambda.1se","Elastic Net","JZS","ZS-null","full model")# ,"ScanBMA")
methods_Interval=c("BMA-bicreg","Spikeslab","UIP", "Benchmark", "BIC ","AIC","EB-local","EB-global",
                   "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","JZS","ZS-null","full model")# ,"ScanBMA")

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


# CRPS.res <- matrix(NA, length(methods_Interval),length(dataname))
# rownames(CRPS.res) <- methods_Interval
# colnames(CRPS.res) <- dataname



for (i in 1:length(filename)){
  load(file=filename[i])
  R2.res[ ,i] <- res1$point.pred.metric[ ,4]
  phat.res[ ,i] <- res1$point.pred.metric[ ,5]
  Coverage.res[ ,i] <- res1$uncertainty.pred.metric[ ,1]
  Width.res[ ,i] <- res1$uncertainty.pred.metric[ ,2]
  MIS.res[ ,i] <- res1$uncertainty.pred.metric[ ,3]
#  CRPS.res[ ,i] <- res1$uncertainty.pred.metric[ ,4]
}



# Divide by JZS prior results 
Coverage.error.res <- abs(Coverage.res -0.95)
Width.res.norm <-  sweep(Width.res, 2, Width.res[ 14,],'/')
MIS.res.norm <-  sweep(MIS.res, 2, MIS.res[14 ,],'/')
#CRPS.res.norm <-  sweep(CRPS.res, 2, CRPS.res[13 ,],'/')





pred_sum_unnorm_alldata <-  list("R2"=R2.res, 
                                 "phat"=phat.res, 
                                 "PredCoverage"=Coverage.res, 
                              "MeanIntScore.predictions"=MIS.res, 
                              #"CRPS.predictions"=CRPS.res,
                              "Width.predictions"=Width.res)

pred_sum_Norm_alldata <-  list("R2"=R2.res, "phat"=phat.res, 
                               "PredCoverageerror"=Coverage.error.res, 
                               "MeanIntScore.pred.norm"=MIS.res.norm, 
                               #"CRPS.pred.norm"=CRPS.res.norm,
                               "Width.pred.norm"=Width.res.norm)


save(pred_sum_unnorm_alldata, file = "./Summary-Rcode/pred_sum_unnorm_alldata.rda")

save(pred_sum_Norm_alldata, file = "./Summary-Rcode/pred_sum_Norm_alldata.rda")



coverage.low <- Coverage.res[ , 1:8]
coverage.high <- Coverage.res[ , 9:12]





# removing EMVS for now because it does not tell me how many variables are selected
# for higher dimensional datasets
R2.res2 <- R2.res
phat.res2 <- phat.res

par(mfrow=c(1,3))
plot(phat.res2[ ,1],R2.res2[ ,1], xlab = "Average model size", ylab = "R2test", main = colnames(R2.res)[1])

plot(phat.res2[ ,2],R2.res2[ ,2], xlab = "Average model size", ylab = "R2test", main = colnames(R2.res)[2])


j=12
cutoff <- min(phat.res2[ ,j]) + 2*sd(phat.res2[ ,j])
par(mfrow=c(1,2))
dat.name <- colnames(R2.res)[j]
gr.title <- paste0(dat.name, " (n=",dat.dim[dat.name,1 ],", p=",dat.dim[dat.name,2 ],")",sep="" )
plot(R2.res2[ ,j],phat.res2[ ,j], ylab = "Average model size", xlab = "Average R2test", 
     main =gr.title,col=ifelse(phat.res2[ ,j]>cutoff,"red","black"),
     cex=2,pch=c(1:18),xlim = c(0,max(R2.res2[ ,j]+0.1)))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

legend(x=-0.5,y=1,legend = rownames(R2.res2),pch=c(1:18),
       col=ifelse(phat.res2[ ,j]>cutoff,"red","black") ,cex=0.8,bty="n",ncol = 2,y.intersp = 3)
# 
# plot(R2.res2[ ,4],phat.res2[ ,4], ylab = "Average model size", xlab = "R2test", 
#   main = colnames(R2.res)[4],col=ifelse(phat.res2[ ,4]>2*min(phat.res2[ ,4]),"red","black"),cex=2,pch=c(1:18))
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("topleft",legend = rownames(R2.res2),pch=c(1:18), ,cex=1,bty="n")

## subset and plot with label
s <- subset(d, In.Deg > 50 & Out.Deg > 50)
with(d, plot(In.Deg, Out.Deg))
text(s[,1:2], labels = s[,3], pos = 1)