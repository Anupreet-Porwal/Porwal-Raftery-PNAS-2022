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
library(SIS)
library(PRROC)

registerDoParallel(cores=2)
#### Functions #### ----


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


cdfBMAbeta <-
  function (x, WEIGHTS, MEAN, SD, offset = 0)
  {
    if(x<0){
      pnorm(x,MEAN,SD)*WEIGHTS-offset
    }
    else if (x>=0){
      pnorm(x,MEAN,SD)*WEIGHTS+(1-WEIGHTS)-offset
    }
    
  }

quantBMAbeta <-
  function(alpha, WEIGHTS, MEAN, SD)
  {
    #
    # copyright 2006-present, University of Washington. All rights reserved.
    # for terms of use, see the LICENSE file
    #
    lower <- MEAN-10*SD
    upper <- MEAN+10*SD
    if(WEIGHTS==0) return(0)
    
    uniroot(cdfBMAbeta, lower = lower, upper = upper,
            WEIGHTS=WEIGHTS, MEAN=MEAN, SD=SD, offset = alpha)$root
  }

meanintscore=function(u,l,x, alpha){
  meanis=0
  for (i in 1:length(x)){
    meanis=meanis+u[i]-l[i]+2/alpha*(l[i]-x[i])*sum(l[i]>x[i])+2/alpha*(x[i]-u[i])*sum(u[i]<x[i])
  }
  meanis=meanis/length(x)
  meanis
}


generate_data=function(fit){
  #Xmat=data[ ,-1]
  
  #nsamp=nrow(Xmat)
  
  #indices <- sample(1:nsamp,size=nsamp,replace=TRUE)
  #Xb= Xmat[indices, ]
  y.hat=fit$fitted.values#[indices]
  #y.true=fit$model[ ,"Y"]
  res.mod=fit$residuals
  Yb= y.hat+ sample(res.mod,length(res.mod), replace = TRUE)
  
  sampy=Yb
}

probne0=function(a){
  sum(abs(a)>0)/length(a)
}

CRPS_VS=function(r, pred.mat, true.mat){
  CRPS=matrix(0, nrow = nrow(pred.mat), 1)
  
  for (v in 1:nrow(pred.mat)){
    for (u in 1:r){
      samp1=sample(pred.mat[v, ], 2, replace = TRUE)
      CRPS[v, 1]=CRPS[v,1]+0.5*abs(samp1[1]-samp1[2])-abs(samp1[1]-true.mat[v])
    }
  }
  CRPS=CRPS/r
  
  return(CRPS)
  
}

auprc <- function(probs, lab){
  if(all(is.na(probs))){
    return (NA)
  }else{
    probs <- probs # exclude the intercept
    fg <- probs[lab==TRUE]
    bg <- probs[lab==FALSE]
    pr <- pr.curve(scores.class0 = fg,scores.class1 = bg)
    
    return(pr$auc.integral)
  }
}

auroc <- function(probs, lab){
  probs <- probs # exclude the intercept
  fg <- probs[lab==TRUE]
  bg <- probs[lab==FALSE]
  roc <- roc.curve(scores.class0 = fg,scores.class1 = bg)
  
  return(roc$auc)
}

conf.scores <- function(x,y,lambda, alpha=1, method){
  ind.mat <- matrix(NA, nrow = length(lambda),ncol = ncol(x))
  colnames(ind.mat) <- colnames(x)
  if(method=="lasso"){
    for (i in 1:length(lambda)){
      mod <- glmnet(x,y,family="gaussian", lambda=lambda[i])
      mod.coef <- coef(mod)
      nzero <- abs(as.matrix(mod.coef)) > 0
      ind.mat[i, ] <- nzero[-1]
    }
    
  }else if(method=="elastic"){
    for (i in 1:length(lambda)){
      mod <- glmnet(x,y,family="gaussian", lambda=lambda[i], alpha=alpha)
      mod.coef <- coef(mod)
      nzero <- abs(as.matrix(mod.coef)) > 0
      ind.mat[i, ] <- nzero[-1]
      
    }
    
    
  }else if(method=="scad"){
    for (i in 1:length(lambda)){
      mod <- ncvreg(x,y,penalty = "SCAD", lambda=lambda[i])
      mod.coef <- coef(mod)
      nzero <- abs(as.matrix(mod.coef)) > 0
      ind.mat[i, ] <- nzero[-1]
      
    }
    
  }else if(method=="mcp"){
    for (i in 1:length(lambda)){
      mod <- ncvreg(x,y,penalty = "MCP", lambda=lambda[i])
      mod.coef <- coef(mod)
      nzero <- abs(as.matrix(mod.coef)) > 0
      ind.mat[i, ] <- nzero[-1]
      
    }
    
  }
  
  return(colMeans(ind.mat))
  
}

Horseshoe.conf.score <- function(betasort){
  alpha <- seq(0.01,0.99,0.01)
  effsamp <- nrow(betasort)
  HS.conf.temp <- matrix(NA, nrow = length(alpha),ncol = ncol(betasort))
  for (i in 1:length(alpha)){
    left <- floor(alpha[i]*effsamp/2)
    right <- ceiling((1-alpha[i]/2)*effsamp)
    left.points <- betasort[left, ]
    right.points <- betasort[right, ]
    HS.conf.temp[i, ] <- as.numeric( 1 - ( (left.points <= 0) & (right.points >= 0) ))
  }
  
  HS.conf.score <- colMeans(HS.conf.temp[rowSums(HS.conf.temp)>0, ])
  
  return(HS.conf.score)
}

sortX<- function(Y, X)
{
  r2vec<- rep(NA, times = ncol(X))
  for (i in 1:ncol(X))
    r2vec[i]<- summary(lm(Y~X[,i]))$r.squared
  initial.order<- order(abs(r2vec),decreasing = TRUE)
  sortedX<- X[, initial.order]
  
  return(list(sortedX = sortedX, initial.order = initial.order))
}

true.model.identification=function(datamat,thresh=0.05){
  
  #  colnames(datamat)[ncol(datamat)]<-"y"
  Y=datamat[ ,ncol(datamat)]
  x=model.matrix(Y~.,data = datamat)[ ,-1]
  
  if(ncol(x)>30){
    sis.mod <- SIS(x,Y,family = "gaussian")
    x.temp <- x[ ,sis.mod$ix]
    if(ncol(x.temp)>30){
      x <- sortX(Y,x.temp)$sortedX[ ,1:30]
    }else {
      x <- sortX(Y,x.temp)$sortedX
    }
    datamat <- cbind(x, Y)
  }
  datamat <- as.data.frame(datamat)
  
  a.subsreg=regsubsets(Y~., data=datamat, nvmax = (ncol(x))
                       , method = "exhaustive", really.big = TRUE)
  a    <- summary(a.subsreg)
  size <- apply(a$which,1,sum)
  names(size) <- NULL
  size <- c(1,size)
  a$which <- a$which[,-1, drop=FALSE]
  a$r2 <- a$rsq
  a$r2 <- pmin(pmax(0, a$r2), 0.999)
  nvar=1
  BIC=matrix(NA,nrow=nrow(a$which),ncol=1)
  for (i in 1:nrow(a$which)){
    x.lm <- cbind.data.frame(Y = Y, as.data.frame(x[, a$which[i,
                                                              , drop = FALSE]]))
    lm.fix <- lm(Y ~ ., data = x.lm)
    BIC[i]=BIC(lm.fix)
    if((!any(summary(lm.fix)$coefficients[2:nrow(summary(lm.fix)$coefficients) ,4]>thresh))& length(lm.fix$coefficients)>=nvar){
      true.mod=lm.fix
      nvar= length(lm.fix$coefficients)
    }
  }
  if(is.null(true.mod)){
    lm.fix <- lm(Y ~ 1, data = datamat)
    true.mod <- lm.fix
  }
  
  return (true.mod)
}

est_varselection_bms=function( datamat, Xmat,true.mod,bootn=100){
  
  tot.num <- ncol(model.matrix(Y~.,datamat))
  truecoef <- matrix(0, nrow = (tot.num), ncol=1)
  rownames(truecoef) <- colnames(model.matrix(Y~.,datamat))
  
  
  for(i in 1:length(truecoef)){
    if(rownames(truecoef)[i] %in% rownames(as.data.frame(true.mod$coefficients))){
      truecoef[i] <- true.mod$coefficients[match(rownames(truecoef)[i],rownames(as.data.frame(true.mod$coefficients)))]
    }
  }
  
  # True coefficients except the intercept
  true.cov.coef=truecoef[2:length(truecoef)]
  
  lvs=c("FALSE","TRUE")
  truecoef.ind <- factor(abs(truecoef)>0, levels = rev(lvs))
  true.cov.coef.ind <- truecoef.ind[2:length(truecoef.ind)]
  
  
  methods_Point=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  methods_Interval=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  
  comp_time_fit=mse_mat_coef=rmse_mat_coef=
    mae_mat_coef=matrix(0, nrow = bootn, ncol = length(methods_Point))
  colnames(comp_time_fit)= colnames(mse_mat_coef)=colnames(mae_mat_coef)=
    colnames(rmse_mat_coef)= methods_Point
  width_mat_coef=coverage_mat_coef=IntervalScore_mat_coef=
    CRPS_mat_coef=array(0, dim=c(bootn,length(methods_Interval),
                                 (length(truecoef)-1))) 
  dimnames(IntervalScore_mat_coef)=dimnames(CRPS_mat_coef)=
    dimnames(coverage_mat_coef)=dimnames(width_mat_coef)=
    list(NULL ,methods_Interval,rownames(truecoef)[2:length(truecoef)])
  
  zero_coverage= nzero_coverage = zero_intscore = nzero_intscore <-
    matrix(0, nrow = bootn, ncol = length(methods_Interval))
  colnames(zero_coverage)=colnames(nzero_coverage)=colnames(zero_intscore)=colnames(nzero_intscore)=methods_Interval
  
  methods_varsel=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  
  sens.mat=spec.mat=prec.mat=F1.mat=type1.mat=
    type2.mat=matrix(0, nrow = bootn, ncol = length(methods_varsel))
  colnames(sens.mat)=colnames(spec.mat)=colnames(type1.mat)=
    colnames(type2.mat)= methods_varsel
  
  methods.auprc=c("EB-local", "g-sqrt n", "Hyper g","JZS-BMA")
  auroc.mat=auprc.mat=matrix(0, nrow = bootn, ncol = length(methods.auprc)) 
  colnames(auroc.mat)=colnames(auprc.mat)=methods.auprc
  
  
  # Y has to be the last column on the dataset! for generate_data function to work.
  
  for (b in 1:bootn){
    print(paste0("Iteration no:",b))
    
    sampdata <- cbind(Xmat, generate_data(true.mod))
    N <- nrow(sampdata)
    colnames(sampdata)[ncol(sampdata)] <- "Y"
    sampdata <- as.data.frame(sampdata)
    
    X.samp <- model.matrix(Y~., sampdata)
    Y.samp <- as.matrix(sampdata[ ,"Y"])
    
    lvs=c("FALSE","TRUE")
    a=0.05
  
    comp_time_fit[b,1] <- system.time({EBlocal.mod <- bas.lm(Y~.,sampdata, prior = "EB-local",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    EBlocal.coef= coef(EBlocal.mod,estimator='HPM')
    EBlocal.conf.coef=confint(EBlocal.coef)[,1:2]
    EBlocal.ind=factor(abs(EBlocal.coef$postmean)>0, levels = rev(lvs))
    
    
    comp_time_fit[b,2] <- system.time({gsqrtn.mod <- bas.lm(Y~., sampdata, prior = "g-prior", alpha = sqrt(N),method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    gsqrtn.coef= coef(gsqrtn.mod,estimator='HPM')
    gsqrtn.conf.coef=confint(gsqrtn.coef)[,1:2]
    gsqrtn.ind=factor(abs(gsqrtn.coef$postmean)>0, levels = rev(lvs))
    
    comp_time_fit[b,3] <- system.time({hypergprior.mod <- bas.lm(Y~.,sampdata, prior = "hyper-g",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    hypergprior.coef= coef(hypergprior.mod,estimator='HPM')
    hypergprior.conf.coef=confint(hypergprior.coef)[,1:2]
    hypergprior.ind=factor(abs(hypergprior.coef$postmean)>0, levels = rev(lvs))
    

    comp_time_fit[b,4] <- system.time({JZS.mod <- bas.lm(Y~.,sampdata, prior = "JZS",method = "MCMC", MCMC.iterations = 10000)})[3]
    JZS.coef= coef(JZS.mod)
    JZS.conf.coef=confint(JZS.coef)[,1:2]
    JZS.probne0=JZS.mod$probne0
    JZS.ind=factor(JZS.probne0>0.5, levels = rev(lvs))
    
    
    coef.matrix=cbind(EBlocal.coef$postmean, 
                      gsqrtn.coef$postmean, 
                      hypergprior.coef$postmean,
                      JZS.coef$postmean)
    coef.matrix <- coef.matrix[-1, ]
    colnames(coef.matrix) <- methods_Point
    
    
    ### MSE MAE for beta's point est
    for (i in 1:ncol(coef.matrix)){
      mse_mat_coef[b,i]=MSE(coef.matrix[ ,i],true.cov.coef)
      rmse_mat_coef[b,i]=rmse(coef.matrix[ ,i],true.cov.coef)
      mae_mat_coef[b,i]=MAE(coef.matrix[ ,i],true.cov.coef)
    }
    
    quantile.list=list( "EBlocal"=EBlocal.conf.coef,
                        "gsqrtn"=gsqrtn.conf.coef,
                        "Hyper g"=hypergprior.conf.coef,
                        "JZS-BMA"=JZS.conf.coef)
    
    
    # Width calculations
    for (j in 1:length(quantile.list)){
      width_mat_coef[b,j, ]= quantile.list[[j]][-1,2]-quantile.list[[j]][-1,1]
    }
    for (e in 1:length(methods_Interval)){
      for (d in 2:length(truecoef)){
        coverage_mat_coef[b,e,d-1]= sum((truecoef[d]>=quantile.list[[e]][d,1])&&(truecoef[d]<=quantile.list[[e]][d,2]))
        IntervalScore_mat_coef[b,e,d-1]=meanintscore(quantile.list[[e]][d,2],quantile.list[[e]][d,1],truecoef[d],a)
      }
    }
    
    zero_coverage[b, ] <- rowMeans(as.matrix(coverage_mat_coef[b,,!as.logical(true.cov.coef.ind)]))
    zero_intscore[b, ] <- rowMeans(as.matrix(IntervalScore_mat_coef[b,,!as.logical(true.cov.coef.ind)]))
    
    nzero_coverage[b, ] <- rowMeans(as.matrix(coverage_mat_coef[b,,as.logical(true.cov.coef.ind)]))
    nzero_intscore[b, ] <- rowMeans(as.matrix(IntervalScore_mat_coef[b,,as.logical(true.cov.coef.ind)]))
    
    
    
    confscore.matrix=data.frame(EBlocal.ind[-1],
                                gsqrtn.ind[-1],
                                hypergprior.ind[-1],
                                JZS.probne0[-1])
    colnames(confscore.matrix)=methods.auprc
    
    for(i in 1:ncol(confscore.matrix)){
      auroc.mat[b,i] <- auroc(confscore.matrix[ , i], true.cov.coef.ind)
      auprc.mat[b,i] <- auprc(confscore.matrix[ , i], true.cov.coef.ind)
    }
    
    ind.matrix=data.frame(EBlocal.ind,
                          gsqrtn.ind, 
                          hypergprior.ind,
                          JZS.ind)
    ind.matrix=ind.matrix[-1, ]
    colnames(ind.matrix)=methods_varsel
    
    
    
    for(i in 1:ncol(ind.matrix)){
      sens.mat[b, i]=sensitivity(ind.matrix[ ,i],true.cov.coef.ind)
      prec.mat[b,i]=posPredValue(ind.matrix[ ,i],true.cov.coef.ind)
      spec.mat[b,i]=specificity(ind.matrix[ ,i],true.cov.coef.ind)
    }
    
    
    type1.mat[b, ]=1-spec.mat[b, ]
    type2.mat[b, ]=1-sens.mat[b, ]
    F1.mat[b, ]=2*prec.mat[b, ]*sens.mat[b, ]/(prec.mat[b, ]+sens.mat[b, ])
    
  }
  
  
  sum_mat=function(mat){
    mat1=rbind(mat, colMean=colMeans(mat))
    mat2=cbind(mat1, rowMean= rowMeans(mat1))
    rownames(mat2)[1:(nrow(mat2)-1)]=rownames(mat)
    colnames(mat2)[1:(ncol(mat2)-1)]=colnames(mat)
    return(mat2)
  }
  
  #### Summary plots and tables ####----
  
  meancoverage_mat_coef=t(apply(coverage_mat_coef, c(2,3), mean))
  sum_meancoverage_mat_coef=sum_mat(meancoverage_mat_coef)
  meanIntervalScore_mat_coef=t(apply(IntervalScore_mat_coef,c(2,3),mean))
  sum_meanIntervalScore_mat_coef=sum_mat(meanIntervalScore_mat_coef)
  meanwidth_mat_coef=t(apply(width_mat_coef,c(2,3),mean))
  sum_meanwidth_mat_coef=sum_mat(meanwidth_mat_coef)
  
  sum_mat_coef=cbind(colMeans(mse_mat_coef), colMeans(mae_mat_coef), colMeans(rmse_mat_coef))
  colnames(sum_mat_coef)=c("AMSE","AMAE","ARMSE")
  sum_mat_coef2=cbind(colMeans(sens.mat),colMeans(spec.mat),colMeans(prec.mat), 
                      colMeans(F1.mat),colMeans(type1.mat), colMeans(type2.mat), 
                      colMeans(type2.mat+type1.mat))
  colnames(sum_mat_coef2)=c("Sensitivity(Recall)", 
                                                    "Specificity",
                                                    "Precision", "F_1 score", 
                                                    "Type I error rate", 
                                                    "Type II error rate", 
                                                    "Total error rate")
  summary_HT=sum_mat(sum_mat_coef2)

  mean_comp_time=as.matrix(t(colMeans(comp_time_fit)))
  
  nr=nrow(sum_meancoverage_mat_coef)
  nc=ncol(sum_meancoverage_mat_coef)
  colsum=rbind(sum_meancoverage_mat_coef[ nr,1:(nc-1)],sum_meanwidth_mat_coef[ nr,1:(nc-1)]
               ,sum_meanIntervalScore_mat_coef[ nr,1:(nc-1)])
  rownames(colsum)=c("Mean coverage","Mean Width","Mean Interval score")
  
  
  results <- list("mse_mat_coef"=mse_mat_coef, "mae_mat_coef"=mae_mat_coef,
                  "rmse_mat_coef"=rmse_mat_coef,
                  "coverage_mat_coef" =apply(coverage_mat_coef,c(1,2),mean),
                  "width_mat_coef"=apply(width_mat_coef,c(1,2),mean),
                  "IntervalScore_mat_coef"=apply(IntervalScore_mat_coef,c(1,2),mean),
                  "sens.mat"=sens.mat,
                  "spec.mat"=spec.mat, "prec.mat"=prec.mat,
                  "F1.mat"=F1.mat, "type1.mat"=type1.mat, "type2.mat"=type2.mat,
                  "total.err.mat"=type1.mat+type2.mat,
                  "pest.metric"=sum_mat_coef, "uncertainty.coef"=colsum,
                  "hypotest"=summary_HT,
                  "comp_time_fit"=comp_time_fit ,"comp.time"=mean_comp_time,
                  "AUPRC.mat"=auprc.mat, "AUROC.mat"=auroc.mat, "zero_coverage"=colMeans(zero_coverage),
                  "nzero_coverage"=colMeans(nzero_coverage), "zero_IntScore"=colMeans(zero_intscore),
                  "nzero_IntScore"=colMeans(nzero_intscore))
  
  return(results)
  
}
