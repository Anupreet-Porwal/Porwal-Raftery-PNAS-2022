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

registerDoParallel(cores=2)
#### Functions #### ----

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

correlation_heatmap=function(cormat){
  upper_tri <- get_upper_tri(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed()+
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  
  print(ggheatmap)
  
}

r2test <- function(a,b){
  SSE <- 0
  SSE <- sum((a-b)^2)
  SST <- sum((a)^2)
  r2 <- 1- SSE/SST
  r2
}


MSE <- function(a,b){
  SSE <- sum((a-b)^2)/length(a)
  SSE
}


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}



MAE <- function(a,b){
  mae<-sum(abs(a-b))/length(a)
  mae
}


cdfBMAnormal <-
  function (x, WEIGHTS, MEAN, SD, offset = 0)
  {
    #
    # copyright 2006-present, University of Washington. All rights reserved.
    # for terms of use, see the LICENSE file
    #
    sum(WEIGHTS*pnorm(x, mean = MEAN, sd = SD)) - offset
  }



quantBMAnormal <-
  function(alpha, WEIGHTS, MEAN, SD)
  {
    #
    # copyright 2006-present, University of Washington. All rights reserved.
    # for terms of use, see the LICENSE file
    #
    lower <- min(MEAN-6*SD)
    upper <- max(MEAN+6*SD)
    
    if (cdfBMAnormal(lower, WEIGHTS, MEAN, SD, 0) > alpha) return(NA)
    if (cdfBMAnormal(upper, WEIGHTS, MEAN, SD, 0) < alpha) return(NA)
    
    uniroot(cdfBMAnormal, lower = lower, upper = upper,
            WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)$root
  }


meanintscore<-function(u,l,x, alpha){
  meanis<-0
  for (i in 1:length(x)){
    meanis<-meanis+u[i]-l[i]+2/alpha*(l[i]-x[i])*sum(l[i]>x[i])+2/alpha*(x[i]-u[i])*sum(u[i]<x[i])
  }
  meanis<-meanis/length(x)
  meanis
}

CRPS_prediction<-function(r, pred.mat, true.mat){
  CRPS<-matrix(0, nrow = nrow(pred.mat), 1)
  
  for (v in 1:nrow(pred.mat)){
    for (u in 1:r){
      samp1<-sample(pred.mat[v, ], 2, replace = TRUE)
      CRPS[v, 1]<-CRPS[v,1]+0.5*abs(samp1[1]-samp1[2])-abs(samp1[1]-true.mat[v])
    }
  }
  CRPS=CRPS/r
  
  return(mean(CRPS))
  
}

post.samples<-function(n.post, bma.obj, test.set ){
  post.samp<-matrix(NA, nrow = nrow(test.set), ncol = n.post)
  
  for (w in 1:ncol(post.samp)){
    index<-sample(1:bma.obj$n.models, 1, replace = TRUE, prob = bma.obj$postprob)
    post.mean<-bma.obj$mle[index, ] %*% t(test.set)
    mod.sd<-sqrt(bma.obj$residvar[index])
    post.samp[ ,w]<-rnorm(nrow(post.samp), mean=post.mean, sd= mod.sd)
  }
  return(post.samp)
}


post.samples.gprior<-function(n.post, gprior.obj,gprior.pred.obj, test.set,estimator='HPM' ){
  post.samp<-matrix(NA, nrow = nrow(test.set), ncol = n.post)
  
  for (w in 1:ncol(post.samp)){
    if(estimator=='HPM'){
      index<-which.max(gprior.obj$postprobs)
      post.mean<-gprior.pred.obj$Ypred
      
    }else if (estimator=='BMA'){
      index<-sample(1:gprior.obj$n.models, 1, replace = TRUE, prob = gprior.obj$postprobs)  
      post.mean<-gprior.pred.obj$Ypred[which(gprior.pred.obj$best==index), ]
    }
    
    mod.sd<-sqrt(gprior.obj$mse[index])
    post.samp[ ,w]<-rnorm(nrow(post.samp), mean=post.mean, sd= mod.sd)
  }
  return(post.samp)
}

residualVariance <- function(coef, x, y)
  sum((y - x %*% coef)^2)/(length(y) - length(coef))



prediction_study_bms=function(Xmat,Y, bootn=100, split=0.75){
  
  #### Loop over Bootstrap samples #### ####----
  tot.num <- ncol(model.matrix(Y~.,cbind(Xmat,Y)))
  
  # methods included in the point/ Interval estimation analysis
  methods_Point=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  methods_Interval=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  
  mse_mat=rmse_mat=mae_mat= r2test.mat=
    matrix(NA,nrow = bootn,ncol = length(methods_Point))
  colnames(mse_mat)=colnames(mae_mat)=colnames(rmse_mat)= colnames(r2test.mat) =methods_Point
  
  width_mat=coverage_mat=IntervalScore_mat=CRPS_mat=
    matrix(NA, nrow = bootn, ncol = length(methods_Interval)) 
  colnames(width_mat)= colnames(coverage_mat)=colnames(IntervalScore_mat)=
    colnames(CRPS_mat)=methods_Interval
  
  phat.mat <- matrix(NA,nrow = bootn,ncol = length(methods_Point))
  colnames(phat.mat) <- methods_Point
  
  
  # For predictive performance, we are not sampling Y,
  #just do a cross validation with the original dataset using various train test splits
  for (b in 1:bootn){
    print(paste0("Iteration no:",b))
    sampdata=cbind(Xmat, Y)
    N=nrow(sampdata)
    colnames(sampdata)[ncol(sampdata)]="Y"
    N.train <- round(N * split)
    N.test <- N - N.train
    
    is.train <- sample(1:N, N.train)
    mean.y.train <- mean(sampdata$Y[is.train])
    sampdata$Y <- sampdata$Y- mean.y.train
    test <- sampdata[-is.train, ]
    train <- sampdata[is.train, ]
    
    X.samp <- model.matrix(Y~., train)
    Y.samp=as.matrix(train[ ,"Y"])
    y.true.test=test[ ,ncol(test)]
    
    test.samp=model.matrix(Y~., test)
    a=0.05
    
    
    #EB-local
    EBlocal.mod <- bas.lm(Y~., train, prior = "EB-local",method = "MCMC",MCMC.iterations = 10000)
    EBlocal.pred <- predict(EBlocal.mod, test[ ,-ncol(test)], se.fit=TRUE,estimator='HPM')
    EBlocal.pred.mat <- post.samples.gprior(2500, EBlocal.mod,EBlocal.pred, test.samp,estimator = 'HPM')
    EBlocal.conf <- confint(EBlocal.pred)
    EBlocal.phat <- EBlocal.mod$size[which.max(EBlocal.mod$postprobs)]
    
    
    #g = sqrt(N)
    gsqrtn.mod <- bas.lm(Y~., train, prior = "g-prior", alpha = sqrt(N.train),method = "MCMC",MCMC.iterations = 10000)
    gsqrtn.pred <- predict(gsqrtn.mod, test[ ,-ncol(test)], se.fit=TRUE,estimator='HPM')
    gsqrtn.pred.mat <- post.samples.gprior(2500, gsqrtn.mod,gsqrtn.pred, test.samp,estimator = 'HPM')
    gsqrtn.conf <- confint(gsqrtn.pred)
    gsqrtn.phat <- gsqrtn.mod$size[which.max(gsqrtn.mod$postprobs)]
    
    
    # Hyper g priors
    hypergprior.mod <- bas.lm(Y~.,train, prior = "hyper-g",method = "MCMC",MCMC.iterations = 10000)
    hypergprior.pred <- predict(hypergprior.mod, test[ ,-ncol(test)], se.fit=TRUE,estimator='HPM')
    hypergprior.pred.mat <- post.samples.gprior(2500, hypergprior.mod, hypergprior.pred, test.samp,estimator = 'HPM')
    hypergprior.conf <- confint(hypergprior.pred)
    hyperg.phat <- hypergprior.mod$size[which.max(hypergprior.mod$postprobs)]
    
    #JZS
    JZS.mod <- bas.lm(Y~., train, prior = "JZS",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
    JZS.pred <- predict(JZS.mod, test[ ,-ncol(test)], se.fit=TRUE)
    JZS.pred.mat <- post.samples.gprior(2500, JZS.mod,JZS.pred, test.samp,estimator = "BMA")
    JZS.conf <- confint(JZS.pred)
    JZS.phat <- sum(JZS.mod$postprobs*JZS.mod$size)
    
    
    pred.mat <- list(EBlocal.pred$fit, 
                     gsqrtn.pred$fit,
                     hypergprior.pred$fit,
                     JZS.pred$fit)
    names(pred.mat) <- methods_Point
    
    for (i in 1:length(methods_Point)){
      mse_mat[b,i] <- MSE(y.true.test, pred.mat[[i]])
      rmse_mat[b,i] <- rmse(y.true.test, pred.mat[[i]])
      mae_mat[b,i] <- MAE(y.true.test, pred.mat[[i]])
      r2test.mat[b,i] <- r2test(y.true.test,pred.mat[[i]])
    }
    
    phat.mat[b, ] <- cbind(EBlocal.phat,gsqrtn.phat, hyperg.phat,JZS.phat)
    
    conf.mat <- list(EBlocal.conf[ ,1:2], 
                     gsqrtn.conf[ ,1:2],
                     hypergprior.conf[ ,1:2],
                     JZS.conf[ ,1:2])
    
    for(j in 1:length(methods_Interval)){
      conf.method <- conf.mat[[j]]
      width_mat[b,j] <- mean(conf.method[ ,2]-conf.method[ ,1])
      coverage_mat[b,j] <- sum(((y.true.test>conf.method[ ,1])&(y.true.test<conf.method[ ,2])))/ length(y.true.test)
      IntervalScore_mat[b,j] <- meanintscore(conf.method[ ,2],conf.method[ ,1], y.true.test, a)
    }
    
    methods.pred.mat <- list(EBlocal.pred.mat,gsqrtn.pred.mat,
                             hypergprior.pred.mat,JZS.pred.mat)
    
    for(j in 1:(length(methods_Interval)-1)){
      CRPS_mat[b,j] <- CRPS_prediction(100, methods.pred.mat[[j]]  ,y.true.test)
      
    }
    
  }
  
  
  
  #### Summary plots and tables #### ####----
  par(mfrow=c(2,4))
  
  for (w in 1:ncol(width_mat)){
    hist(width_mat[ ,w], main = paste("PI-",colnames(width_mat)[w]), xlab = "Avg Width")
  }
  
  
  sum_mat=cbind(colMeans(mse_mat), colMeans(mae_mat), colMeans(rmse_mat),
                colMeans(r2test.mat),colMeans(phat.mat))
  sum_mat2=cbind(colMeans(coverage_mat), colMeans(width_mat), colMeans(IntervalScore_mat),
                 colMeans(CRPS_mat))
  colnames(sum_mat)=c("AMSE","AMAE","ARMSE","Average R^2 test", "average model size")
  colnames(sum_mat2)=c("Mean coverage","Mean Width","Mean Interval score","CRPS")
  
  
  results <- list("mse_mat"=mse_mat,"mae_mat"=mae_mat,"rmse_mat"=rmse_mat, "r2test_mat" = r2test.mat,
                  "phat.mat"=phat.mat,"coverage_mat"=coverage_mat, "width_mat"=width_mat, "CRPS_mat"=CRPS_mat,
                  "IntervalScore"=IntervalScore_mat, "point.pred.metric"=sum_mat,
                  "uncertainty.pred.metric"=sum_mat2)
  return (results)
}
