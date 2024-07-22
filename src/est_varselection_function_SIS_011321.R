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

est_varselection_study=function( datamat, Xmat,true.mod,bootn=100){

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


  methods_Point=c("BMA-bicreg","Spikeslab", "UIP", "BIC","AIC","EB-local","EB-global",
                  "g-sqrt n", "g-1","Hyper g","NLP","Horseshoe","SS Lasso","EMVS"
                  ,"SCAD","MCP","LASSO-lambda.min","LASSO-lambda.1se","Elastic net",
                  "JZS","ZS-null" , "true.mod", "full model")
  methods_Interval=c("BMA -bicreg","Spikeslab","UIP", "BIC", "AIC","EB-local",
                     "EB-global","g-sqrt n" , "g-1","Hyper g", "NLP", "Horseshoe",
                     "JZS","ZS-null" , "true.mod", "full model")#,"ScanBMA")

  comp_time_fit=mse_mat_coef=rmse_mat_coef=
    mae_mat_coef=matrix(0, nrow = bootn, ncol = length(methods_Point))
  colnames(comp_time_fit)= colnames(mse_mat_coef)=colnames(mae_mat_coef)=
    colnames(rmse_mat_coef)= methods_Point
  width_mat_coef=coverage_mat_coef=IntervalScore_mat_coef=
    CRPS_mat_coef=array(0, dim=c(bootn,length(methods_Interval),
                                 (length(truecoef)-1))) # EMVS SS LASSO, LASSO, SCAD, MCP excluded
    dimnames(IntervalScore_mat_coef)=dimnames(CRPS_mat_coef)=
      dimnames(coverage_mat_coef)=dimnames(width_mat_coef)=
      list(NULL ,methods_Interval,rownames(truecoef)[2:length(truecoef)])

    zero_coverage= nzero_coverage = zero_intscore = nzero_intscore <-
      matrix(0, nrow = bootn, ncol = length(methods_Interval))
    colnames(zero_coverage)=colnames(nzero_coverage)=colnames(zero_intscore)=colnames(nzero_intscore)=methods_Interval

    methods_varsel=c("BMA-bicreg","Spikeslab", "UIP", "BIC","AIC","EB-local",
                     "EB-global","g-sqrt n", "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS",
                     "SCAD","MCP","LASSO-lambda.min","LASSO-lambda.1se","Elastic net",
                     "JZS","ZS-null", "true.mod" , "full model")
    sens.mat5=spec.mat5=prec.mat5=F1.mat5=type1.mat5=
    type2.mat5=matrix(0, nrow = bootn, ncol = length(methods_varsel))
  colnames(sens.mat5)=colnames(spec.mat5)=colnames(type1.mat5)=
    colnames(type2.mat5)= methods_varsel

  sens.mat9=spec.mat9=prec.mat9=F1.mat9=type1.mat9=type2.mat9=
    matrix(0, nrow = bootn, ncol = length(methods_varsel))
  colnames(sens.mat9)=colnames(spec.mat9)=colnames(type1.mat9)=colnames(type2.mat9)= methods_varsel

  methods.auprc=c("BMA-bicreg","Spikeslab", "UIP", "BIC","AIC","EB-local",
                  "EB-global","g-sqrt n", "g-1","Hyper g","NLP", "Horseshoe", "SS Lasso","EMVS",
                  "SCAD","MCP","LASSO","Elastic net", "JZS","ZS-null", "true.mod" , "full model")
  auroc.mat=auprc.mat=matrix(0, nrow = bootn, ncol = length(methods.auprc)) # removing second lasso
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

    if(ncol(X.samp)>30){
      comp.time.sis <- system.time({sis.mod <- SIS(X.samp[ ,-1],Y.samp,family = "gaussian")})[3]
      X.samp.scr <- cbind(rep(1,N), X.samp[ ,sis.mod$ix+1])
      colnames(X.samp.scr) <- c("(Intercept)",colnames(X.samp)[sis.mod$ix+1])
      sampdata.scr <- cbind(X.samp.scr[ ,-1], Y.samp)
      colnames(sampdata.scr) <- c(colnames(X.samp)[sis.mod$ix+1],"Y")

      # Selecting top 100 based on R^2
      # X.samp.scr <-  cbind(rep(1,N),sortX(Y.samp,X.samp[ ,-1])$sortedX[ ,1:100])
      # colnames(X.samp.scr)[1] <- "(Intercept)"
      # sampdata.scr <- cbind(X.samp.scr[ ,-1], Y.samp)
      # colnames(sampdata.scr)[ncol(sampdata.scr)] <- "Y"
    }
    sampdata.scr <- as.data.frame(sampdata.scr)


    comp_time_fit[b,1] <- system.time({
      if (ncol(X.samp.scr)>30){
      ibma.mod <- iBMA.bicreg(subset(sampdata.scr,select=-Y), sampdata.scr[ ,ncol(sampdata.scr)], strict = FALSE,OR=20)
      bma.mod <- ibma.mod$bma
      } else{
        bma.mod <- bicreg(subset(sampdata.scr,select=-Y), sampdata.scr[ ,ncol(sampdata.scr)], strict = FALSE,OR=100)
      }
      })[3]+comp.time.sis
    bma.coef.scr=bma.mod$postmean
    a=0.05
    bma.probne0.scr=c(100, bma.mod$probne0)
    #scanbma.probne0.scr=c(100, scanbma.mod$probne0)
    bma.condpostmean.scr <- bma.mod$condpostmean
    bma.condpostsd.scr <- bma.mod$condpostsd
    #scanbma.condpostmean.scr <- scanbma.mod$condpostmean
    #scanbma.condpostsd.scr <- scanbma.mod$condpostsd
    #names(scanbma.probne0.scr) = names(scanbma.condpostmean.scr) = names( scanbma.condpostsd.scr)
    names(bma.probne0.scr) =names(bma.condpostmean.scr)=
      names(bma.condpostsd.scr) <- c("(Intercept)",bma.mod$input.names)

    #scanbma.coef= scanbma.probne0 = scanbma.condpostmean = scanbma.condpostsd =
    bma.coef = bma.probne0 = bma.condpostmean =
      bma.condpostsd  <- matrix(0, nrow = nrow(truecoef), ncol = 1)
    #rownames(scanbma.coef)= rownames(scanbma.probne0) = rownames(scanbma.condpostmean)
    #= rownames(scanbma.condpostsd)=
    rownames(bma.coef)=rownames(bma.probne0)=
      rownames(bma.condpostmean)=rownames(bma.condpostsd) <- rownames(truecoef)


    for(q in 1:nrow(bma.coef)){
      if(rownames(bma.coef)[q] %in% names(bma.coef.scr)){
        bma.coef[q,1]=bma.coef.scr[match(rownames(bma.coef)[q],names(bma.coef.scr))]
        #scanbma.coef[q,1]=scanbma.coef.scr[match(rownames(scanbma.coef)[q],names(scanbma.coef.scr))]
        bma.probne0[q,1]=bma.probne0.scr[match(rownames(bma.probne0)[q],names(bma.probne0.scr))]
        #scanbma.probne0[q,1]=scanbma.probne0.scr[match(rownames(scanbma.probne0)[q],names(scanbma.probne0.scr))]
        bma.condpostmean[q,1]=bma.condpostmean.scr[match(rownames(bma.condpostmean)[q],names(bma.condpostmean.scr))]
        #scanbma.condpostmean[q,1]=scanbma.condpostmean.scr[match(rownames(scanbma.condpostmean)[q],names(scanbma.condpostmean.scr))]
        bma.condpostsd[q,1]=bma.condpostsd.scr[match(rownames(bma.condpostsd)[q],names(bma.condpostsd.scr))]
        #scanbma.condpostsd[q,1]=scanbma.condpostsd.scr[match(rownames(scanbma.condpostsd)[q],names(scanbma.condpostsd.scr))]
      }
    }



    bma.conf.coef=matrix(0,nrow = length(bma.coef), ncol=2)
    rownames(bma.conf.coef) <- rownames(truecoef)
    for (l in 1:length(bma.coef)){

      low=quantBMAbeta(a/2, bma.probne0[l]/100, bma.condpostmean[l],bma.condpostsd[l])
      up=quantBMAbeta(1-a/2, bma.probne0[l]/100, bma.condpostmean[l],bma.condpostsd[l])
      bma.conf.coef[l, ]=c(low,up)
    }
    # 95% CI for each of the coefficients
    # scanbma.conf.coef=matrix(0,nrow = length(scanbma.coef), ncol=2)
    # rownames(scanbma.conf.coef) <- rownames(truecoef)
    # for (l in 1:length(scanbma.coef)){
    #   low=quantBMAbeta(a/2, scanbma.probne0[l]/100, scanbma.condpostmean[l],scanbma.condpostsd[l])
    #   up=quantBMAbeta(1-a/2, scanbma.probne0[l]/100, scanbma.condpostmean[l],scanbma.condpostsd[l])
    #   scanbma.conf.coef[l, ]=c(low,up)
    # }
    #
    bma.ind5=factor(bma.probne0/100> 0.5, levels = rev(lvs))
    bma.ind9=factor(bma.probne0/100> 0.95, levels = rev(lvs))




    comp_time_fit[b,2] <- system.time({ss.mod <- quiet(lm.spike(Y~.,niter = 4000, data = sampdata, error.distribution = "gaussian"))})[3]
    ss.betasamp=ss.mod$beta[1000:nrow(ss.mod$beta), ]
    ss.coef=colMeans(ss.betasamp)
    ss.conf.coef=t(apply(ss.betasamp, 2,quantile, probs=c(0.025,0.975)))
    ss.probne0=apply(ss.betasamp,2,probne0)
    ss.ind5= factor(ss.probne0>0.5, levels = rev(lvs))
    ss.ind9= factor(ss.probne0>0.95, levels = rev(lvs))

    comp_time_fit[b,3] <- system.time({gpriorUnit.mod <- bas.lm(Y~., sampdata, prior = "g-prior", alpha = N,method = "MCMC", MCMC.iterations = 10000)})[3]
    gpriorUnit.coef=coef(gpriorUnit.mod)
    gpriorUnit.conf.coef=confint(gpriorUnit.coef)[ ,1:2]
    gpriorUnit.probne0=gpriorUnit.mod$probne0
    gpriorUnit.ind5=factor(gpriorUnit.probne0>0.5, levels = rev(lvs))
    gpriorUnit.ind9=factor(gpriorUnit.probne0>0.95, levels = rev(lvs))


    comp_time_fit[b,4] <- system.time({gpriorBIC.mod <- bas.lm(Y~., sampdata, prior = "BIC",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    gpriorBIC.coef= coef(gpriorBIC.mod)
    gpriorBIC.conf.coef=confint(gpriorBIC.coef)[,1:2]
    gpriorBIC.probne0=gpriorBIC.mod$probne0
    gpriorBIC.ind5=factor(gpriorBIC.probne0>0.5, levels = rev(lvs))
    gpriorBIC.ind9=factor(gpriorBIC.probne0>0.95, levels = rev(lvs))


    comp_time_fit[b,5] <- system.time({AIC.mod <- bas.lm(Y~.,sampdata, prior = "AIC",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    AIC.coef= coef(AIC.mod)
    AIC.conf.coef=confint(AIC.coef)[,1:2]
    AIC.probne0=AIC.mod$probne0
    AIC.ind5=factor(AIC.probne0>0.5, levels = rev(lvs))
    AIC.ind9=factor(AIC.probne0>0.95, levels = rev(lvs))



    comp_time_fit[b,6] <- system.time({EBlocal.mod <- bas.lm(Y~.,sampdata, prior = "EB-local",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    EBlocal.coef= coef(EBlocal.mod)
    EBlocal.conf.coef=confint(EBlocal.coef)[,1:2]
    EBlocal.probne0=EBlocal.mod$probne0
    EBlocal.ind5=factor(EBlocal.probne0>0.5, levels = rev(lvs))
    EBlocal.ind9=factor(EBlocal.probne0>0.95, levels = rev(lvs))




    comp_time_fit[b,7] <- system.time({EBglobal.mod <- bas.lm(Y~.,sampdata, prior = "EB-global",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    EBglobal.coef= coef(EBglobal.mod)
    EBglobal.conf.coef=confint(EBglobal.coef)[,1:2]
    EBglobal.probne0=EBglobal.mod$probne0
    EBglobal.ind5=factor(EBglobal.probne0>0.5, levels = rev(lvs))
    EBglobal.ind9=factor(EBglobal.probne0>0.95, levels = rev(lvs))


    comp_time_fit[b,8] <- system.time({gsqrtn.mod <- bas.lm(Y~., sampdata, prior = "g-prior", alpha = sqrt(N),method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    gsqrtn.coef= coef(gsqrtn.mod)
    gsqrtn.conf.coef=confint(gsqrtn.coef)[,1:2]
    gsqrtn.probne0=gsqrtn.mod$probne0
    gsqrtn.ind5=factor(gsqrtn.probne0>0.5, levels = rev(lvs))
    gsqrtn.ind9=factor(gsqrtn.probne0>0.95, levels = rev(lvs))

    comp_time_fit[b,9] <- system.time({g1.mod <- bas.lm(Y~., sampdata, prior = "g-prior", alpha = 1,method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    g1.coef= coef(g1.mod)
    g1.conf.coef=confint(g1.coef)[,1:2]
    g1.probne0=g1.mod$probne0
    g1.ind5=factor(g1.probne0>0.5, levels = rev(lvs))
    g1.ind9=factor(g1.probne0>0.95, levels = rev(lvs))

    comp_time_fit[b,10] <- system.time({hypergprior.mod <- bas.lm(Y~.,sampdata, prior = "hyper-g",method = "MCMC", MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N-2))})[3]
    hypergprior.coef= coef(hypergprior.mod)
    hypergprior.conf.coef=confint(hypergprior.coef)[,1:2]
    hypergprior.probne0=hypergprior.mod$probne0
    hypergprior.ind5=factor(hypergprior.probne0>0.5, levels = rev(lvs))
    hypergprior.ind9=factor(hypergprior.probne0>0.95, levels = rev(lvs))


    comp_time_fit[b,11] <- system.time({nlp.modselect=quiet(modelSelection(y=Y.samp,x=X.samp, data= sampdata, scale = FALSE,center = FALSE))})[3]
    nlp.postsamp=rnlp(msfit=nlp.modselect, niter = 3000,burnin = 1000)
    nlp.betasamp=nlp.postsamp[ , 2:(ncol(nlp.postsamp)-1)]
    nlp.coef= colMeans(nlp.betasamp)
    nlp.conf.coef=t(apply(nlp.betasamp, 2,quantile, probs=c(0.025,0.975)))
    NLP.probne0=nlp.modselect$margpp
    NLP.ind5=factor(NLP.probne0>0.5, levels = rev(lvs))
    NLP.ind9=factor(NLP.probne0>0.95, levels = rev(lvs))

    comp_time_fit[b,12] <- system.time({
      HS.mod <- quiet(horseshoe(Y.samp,X.samp,method.tau = "halfCauchy", method.sigma = "Jeffreys",alpha = 0.05))
      alpha <- 0.5
      effsamp <- ncol(HS.mod$BetaSamples)
      left.50 <- floor(alpha*effsamp/2)
      right.50 <- ceiling((1-alpha/2)*effsamp)
      BetaSort <- apply(HS.mod$BetaSamples, 1, sort, decreasing = F)
      left.points <- BetaSort[left.50, ]
      right.points <- BetaSort[right.50, ]
      HS.vs5 <- as.numeric( 1 - ( (left.points <= 0) & (right.points >= 0) ))
      HS.vs9 <- HS.var.select(HS.mod, Y.samp, method = "intervals")
    })[3]
    HS.coef=HS.mod$BetaHat
    HS.conf.coef=cbind(HS.mod$LeftCI, HS.mod$RightCI)
    HS.conf.score <- Horseshoe.conf.score(BetaSort[ ,-1])
    HS.ind5=factor(HS.vs5)
    HS.ind9=factor(HS.vs9)
    HS.ind5= revalue(HS.ind5, c("0"="FALSE", "1"="TRUE"))
    HS.ind9= revalue(HS.ind9, c("0"="FALSE", "1"="TRUE"))

    comp_time_fit[b,13] <- system.time({sslasso.mod <- SSLASSO(X.samp[ ,-1], Y.samp, variance = "unknown")})[3] # you dont need to give a 1s column seperately
    sslasso.coef=cbind(sslasso.mod$intercept[100],t(sslasso.mod$beta[ ,100 ]))
    sslasso.conf.score <- rowMeans(sslasso.mod$select)
    sslasso.ind=factor(abs(sslasso.coef)>0, levels = rev(lvs))


    #EMVS
    comp_time_fit[b,14] <- system.time({
      sampdata.EM=X.samp[ ,-1]#matrix(as.numeric(unlist(sampdata)),nrow=nrow(sampdata),ncol=ncol(sampdata))
      sampdata.EM=cbind(sampdata.EM,Y.samp)
      colnames(sampdata.EM)=c(colnames(X.samp)[-1],"Y")
      p=ncol(sampdata.EM)-1
      epsilon=10^{-5}
      a.em=b.em=1
      beta.init=rep(1,p)
      sigma.init=1
      v0=seq(0.001,10,length.out=1000)
      v1=100
      # independent =TRUE, temp=1 and betainit set to default
      EMVS.mod=quiet(EMVS(Y= sampdata.EM[ ,ncol(sampdata.EM)], X = sampdata.EM[ ,-ncol(sampdata.EM)],v0=v0,v1=v1,
                          type="betabinomial",
                          sigma_init = sigma.init ,
                          epsilon=epsilon,a=a.em,b=b.em))
      # EMVS.mod=quiet(EMVS(Y= sampdata.EM[ ,ncol(sampdata.EM)], X = sampdata.EM[ ,-ncol(sampdata.EM)],v0=v0,v1=v1,
      #                     type="betabinomial",beta_init = beta.init,
      #                     sigma_init = sigma.init ,independent = FALSE,
      #                     epsilon=epsilon,a=a.em,b=b.em, temperature=0.1))
      #EMVS.best=quiet(EMVSbest(EMVS.mod))
    })[3]
    #comp_time_fit[b,15] <- system.time({scanbma.mod=ScanBMA(as.matrix(X.samp.scr[ ,-1]), Y.samp)})[3]+comp.time.sis
    emvs.coef=c(0,EMVS.mod$betas[1, ])
    emvs.probne0=c(1,EMVS.mod$prob_inclusion[1, ])
    emvs.ind5=factor(emvs.probne0>0.5, levels = rev(lvs))
    emvs.ind9=factor(emvs.probne0>0.95, levels = rev(lvs))


    comp_time_fit[b,15] <- system.time({ scad.mod <-  cv.ncvreg(X.samp[ ,-1], Y.samp, family = "gaussian", penalty = "SCAD")})[3]
    #scanbma.coef.scr=scanbma.mod$postmean
    #names(scanbma.coef.scr) <- colnames(X.samp.scr)
    scad.coef=coef(scad.mod)
    scad.conf.score <- conf.scores(X.samp[ ,-1], Y.samp, lambda = scad.mod$lambda, method = "scad")
    scad.ind=factor(abs(as.matrix(scad.coef))>0, levels = rev(lvs))



    comp_time_fit[b,16] <- system.time({ mcp.mod <-  cv.ncvreg(X.samp[ ,-1], Y.samp, family = "gaussian", penalty = "MCP")})[3]
    mcp.coef=coef(mcp.mod)
    mcp.conf.score <- conf.scores(X.samp[ ,-1], Y.samp, lambda = mcp.mod$lambda, method = "mcp")
    mcp.ind=factor(abs(as.matrix(mcp.coef))>0, levels = rev(lvs))


    comp_time_fit[b,17] <- system.time({lasso.mod <- cv.glmnet(X.samp[ , -1], Y.samp, family = "gaussian")})[3]
    lasso.coef=coef(lasso.mod, s="lambda.1se")
    lasso.coef.min=coef(lasso.mod, s="lambda.min")
    lasso.conf.score <- conf.scores(X.samp[ ,-1], Y.samp, lambda = lasso.mod$lambda, method = "lasso")
    lasso.ind=factor(abs(as.matrix(lasso.coef))>0, levels = rev(lvs))
    lasso.ind.min=factor(abs(as.matrix(lasso.coef.min))>0, levels = rev(lvs))


    comp_time_fit[b,18] <- comp_time_fit[b,17]



    comp_time_fit[b,19] <- system.time({
      alp <- seq(0, 0.99, 0.04)
      search <- foreach(i = alp, .combine = rbind,.packages = "glmnet") %dopar% {
        cv <- cv.glmnet(X.samp[ , -1], Y.samp, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
        list(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i, lambda=t(cv$lambda))
      }
      cvm=unlist(search[c(1:25)])
      en.lambda.min <- unlist(search[c(26:50)])
      en.alpha <- unlist(search[c(51:75)])
      min.ind <- which.min(cvm)
      #cv3 <- search[search$cvm == min(search$cvm), ]
      elastic.mod <- glmnet(X.samp[ , -1], Y.samp, family = "gaussian", lambda = en.lambda.min[min.ind], alpha = en.alpha[min.ind])
    })[3]
    elastic.coef <- coef(elastic.mod)
    elastic.conf.score <- conf.scores(X.samp[ ,-1], Y.samp, lambda = search[[75+min.ind]], alpha = en.alpha[min.ind], method = "elastic")
    elastic.ind=factor(abs(as.matrix(elastic.coef))>0, levels = rev(lvs))


    #scanbma.ind5=factor(scanbma.probne0/100> 0.5, levels = rev(lvs))
    #scanbma.ind9=factor(scanbma.probne0/100> 0.95, levels = rev(lvs))
    comp_time_fit[b,20] <- system.time({JZS.mod <- bas.lm(Y~.,sampdata, prior = "JZS",method = "MCMC", MCMC.iterations = 10000)})[3]
    JZS.coef= coef(JZS.mod)
    JZS.conf.coef=confint(JZS.coef)[,1:2]
    JZS.probne0=JZS.mod$probne0
    JZS.ind5=factor(JZS.probne0>0.5, levels = rev(lvs))
    JZS.ind9=factor(JZS.probne0>0.95, levels = rev(lvs))

    comp_time_fit[b,21] <- system.time({ZSnull.mod <- bas.lm(Y~.,sampdata, prior = "ZS-null",method = "MCMC", MCMC.iterations = 10000)})[3]
    ZSnull.coef= coef(ZSnull.mod)
    ZSnull.conf.coef=confint(ZSnull.coef)[,1:2]
    ZSnull.probne0=ZSnull.mod$probne0
    ZSnull.ind5=factor(ZSnull.probne0>0.5, levels = rev(lvs))
    ZSnull.ind9=factor(ZSnull.probne0>0.95, levels = rev(lvs))


    X.true <- X.samp[ ,-1]
    X.true <- X.true[ , as.logical(true.cov.coef.ind)]
    datamat.true <- cbind(X.true,Y.samp)
    colnames(datamat.true)[ ncol(datamat.true)] <- "Y"
    datamat.true <- as.data.frame(datamat.true)
    comp_time_fit[b,22] <- system.time({t.mod <- lm(Y~., data = datamat.true)})[3]
    tmod.coef.scr <- coef(t.mod)
    tmod.conf.coef.scr <- confint(t.mod)
    tmod.coef <- matrix(0, nrow=p+1,ncol = 1)
    tmod.probne0 <- matrix(0, nrow=p+1,ncol = 1)
    tmod.conf.coef <- matrix(0, nrow=p+1,ncol = 2)
    colnames(tmod.conf.coef) <- colnames(tmod.conf.coef.scr)
    cz=1
    for(z in 1:ncol(X.samp)){
      if(truecoef.ind[z]==TRUE){
        tmod.coef[z]=tmod.coef.scr[cz]
        tmod.conf.coef[z, ]=tmod.conf.coef.scr[cz,]
        tmod.probne0[z]=1
        cz=cz+1
      }else if(truecoef.ind[z]==FALSE){next}
    }
    tmod.ind5=factor(tmod.probne0>0.5, levels = rev(lvs))
    tmod.ind9=factor(tmod.probne0>0.95, levels = rev(lvs))

    if(p<N){
    comp_time_fit[b,23] <- system.time({full.mod <- lm(Y~.,data = sampdata)})[3]
    fullmod.coef <- coef(full.mod)
    fullmod.conf.coef <- confint(full.mod)
    fullmod.probne0 <- rep(1,p+1)
    fullmod.ind5=factor(fullmod.probne0>0.5, levels = rev(lvs))
    fullmod.ind9=factor(fullmod.probne0>0.95, levels = rev(lvs))
    }else{
      fullmod.coef <- rep(NA,p+1)
      fullmod.conf.coef <-matrix(NA,nrow=p+1,ncol=2)
      fullmod.probne0 <- rep(NA,p+1)
      fullmod.ind5=factor(fullmod.probne0>0.5, levels = rev(lvs))
      fullmod.ind9=factor(fullmod.probne0>0.95, levels = rev(lvs))
    }



    coef.matrix=cbind(bma.coef[-1],ss.coef[-1],gpriorUnit.coef$postmean[-1],gpriorBIC.coef$postmean[-1],
                      AIC.coef$postmean[-1],EBlocal.coef$postmean[-1], EBglobal.coef$postmean[-1],
                      gsqrtn.coef$postmean[-1], g1.coef$postmean[-1],hypergprior.coef$postmean[-1],
                      nlp.coef[-1],HS.coef[-1],sslasso.coef[-1],emvs.coef[-1],#scanbma.coef[-1],
                      scad.coef[-1],mcp.coef[-1],lasso.coef.min[-1],lasso.coef[-1],elastic.coef[-1],
                      JZS.coef$postmean[-1],ZSnull.coef$postmean[-1],tmod.coef[-1] ,fullmod.coef[-1])
    colnames(coef.matrix) <- methods_Point

    ### MSE MAE for beta's point est
    for (i in 1:ncol(coef.matrix)){
      mse_mat_coef[b,i]=MSE(coef.matrix[ ,i],true.cov.coef)
      rmse_mat_coef[b,i]=rmse(coef.matrix[ ,i],true.cov.coef)
      mae_mat_coef[b,i]=MAE(coef.matrix[ ,i],true.cov.coef)
    }

    quantile.list=list("BMA"=bma.conf.coef,"SpikeSlab"=ss.conf.coef,"UIP"=gpriorUnit.conf.coef,
                       "BIC"=gpriorBIC.conf.coef,"AIC"=AIC.conf.coef, "EBlocal"=EBlocal.conf.coef,
                       "EBglobal"=EBglobal.conf.coef,"gsqrtn"=gsqrtn.conf.coef,"g1"=g1.conf.coef ,
                       "Hyper g"=hypergprior.conf.coef,"NLP"=nlp.conf.coef,"HS"=HS.conf.coef
                       ,"JZS"=JZS.conf.coef, "ZSnull"=ZSnull.conf.coef,"truemod"=tmod.conf.coef
                       ,"fullmod"=fullmod.conf.coef)#,"scanBMA"=scanbma.conf.coef)


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



    confscore.matrix=data.frame(bma.probne0[-1]/100,ss.probne0[-1],gpriorUnit.probne0[-1],gpriorBIC.probne0[-1],
                            AIC.probne0[-1],EBlocal.probne0[-1],EBglobal.probne0[-1],
                           gsqrtn.probne0[-1],g1.probne0[-1],hypergprior.probne0[-1], NLP.probne0[-1],
                           HS.conf.score,sslasso.conf.score,emvs.probne0[-1], #scanbma.ind5,
                           scad.conf.score,mcp.conf.score,lasso.conf.score,elastic.conf.score,
                           JZS.probne0[-1],ZSnull.probne0[-1],tmod.probne0[-1],fullmod.probne0[-1])
    #confscore.matrix=confscore.matrix[-1, ]
    colnames(confscore.matrix)=methods.auprc

    for(i in 1:ncol(confscore.matrix)){
      auroc.mat[b,i] <- auroc(confscore.matrix[ , i], true.cov.coef.ind)
      auprc.mat[b,i] <- auprc(confscore.matrix[ , i], true.cov.coef.ind)
    }

    ind5.matrix=data.frame(bma.ind5,ss.ind5,gpriorUnit.ind5,gpriorBIC.ind5,
                          AIC.ind5,EBlocal.ind5,EBglobal.ind5,gsqrtn.ind5,g1.ind5, hypergprior.ind5,
                           NLP.ind5,HS.ind9,sslasso.ind, emvs.ind5, #scanbma.ind5, # Horseshoe 95% Credible interval are used as a similar analogue to 0.5 threshold of PIP
                           scad.ind,mcp.ind,lasso.ind.min,lasso.ind,elastic.ind, JZS.ind5,ZSnull.ind5,
                           tmod.ind5,fullmod.ind5)
    ind5.matrix=ind5.matrix[-1, ]
    colnames(ind5.matrix)=methods_varsel

    ind9.matrix=data.frame(bma.ind9,ss.ind9,gpriorUnit.ind9,gpriorBIC.ind9,
                           AIC.ind9,EBlocal.ind9,EBglobal.ind9, gsqrtn.ind9,g1.ind5,hypergprior.ind9,
                           NLP.ind9,HS.ind5,sslasso.ind,emvs.ind9,#scanbma.ind9,# Horseshoe 50% Credible interval are used as a similar analogue to 0.95 threshold of PIP
                           scad.ind,mcp.ind,lasso.ind.min,lasso.ind,elastic.ind, JZS.ind9,ZSnull.ind9,
                           tmod.ind9,fullmod.ind9)
    ind9.matrix=ind9.matrix[-1, ]
    colnames(ind9.matrix)=methods_varsel


    for(i in 1:ncol(ind5.matrix)){
      sens.mat5[b, i]=sensitivity(ind5.matrix[ ,i],true.cov.coef.ind)
      sens.mat9[b,i]=sensitivity(ind9.matrix[ ,i],true.cov.coef.ind)
      prec.mat5[b,i]=posPredValue(ind5.matrix[ ,i],true.cov.coef.ind)
      prec.mat9[b,i]=posPredValue(ind9.matrix[ ,i],true.cov.coef.ind)
      spec.mat5[b,i]=specificity(ind5.matrix[ ,i],true.cov.coef.ind)
      spec.mat9[b,i]=specificity(ind9.matrix[ ,i],true.cov.coef.ind)
    }


    type1.mat5[b, ]=1-spec.mat5[b, ]
    type2.mat5[b, ]=1-sens.mat5[b, ]
    type1.mat9[b, ]=1-spec.mat9[b, ]
    type2.mat9[b, ]=1-sens.mat9[b, ]
    F1.mat5[b, ]=2*prec.mat5[b, ]*sens.mat5[b, ]/(prec.mat5[b, ]+sens.mat5[b, ])
    F1.mat9[b, ]=2*prec.mat9[b, ]*sens.mat9[b, ]/(prec.mat9[b, ]+sens.mat9[b, ])

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
  sum_mat_coef2=cbind(colMeans(sens.mat5),colMeans(spec.mat5),colMeans(prec.mat5), colMeans(F1.mat5),colMeans(type1.mat5), colMeans(type2.mat5), colMeans(type2.mat5+type1.mat5))
  sum_mat_coef3=cbind(colMeans(sens.mat9),colMeans(spec.mat9),colMeans(prec.mat9), colMeans(F1.mat9), colMeans(type1.mat9), colMeans(type2.mat9), colMeans(type2.mat9+type1.mat9))
  colnames(sum_mat_coef3)=colnames(sum_mat_coef2)=c("Sensitivity(Recall)", "Specificity","Precision", "F_1 score", "Type I error rate", "Type II error rate", "Total error rate")
  summary_HT_05=sum_mat(sum_mat_coef2)
  summary_HT_95=sum_mat(sum_mat_coef3)

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
                  "sens.mat9"=sens.mat9, "spec.mat9"=spec.mat9, "prec.mat9"=prec.mat9,
                  "F1.mat9"=F1.mat9, "type1.mat9"=type1.mat9, "type2.mat9"=type2.mat9,
                  "total.err.mat9"=type1.mat9+type2.mat9, "sens.mat5"=sens.mat5,
                  "spec.mat5"=spec.mat5, "prec.mat5"=prec.mat5,
                  "F1.mat5"=F1.mat5, "type1.mat5"=type1.mat5, "type2.mat5"=type2.mat5,
                  "total.err.mat5"=type1.mat5+type2.mat5,
                  "pest.metric"=sum_mat_coef, "uncertainty.coef"=colsum,
                  "hypotest.95"=summary_HT_95,"hypotest.50"=summary_HT_05,
                  "comp_time_fit"=comp_time_fit ,"comp.time"=mean_comp_time,
                  "AUPRC.mat"=auprc.mat, "AUROC.mat"=auroc.mat, "zero_coverage"=colMeans(zero_coverage),
                  "nzero_coverage"=colMeans(nzero_coverage), "zero_IntScore"=colMeans(zero_intscore),
                  "nzero_IntScore"=colMeans(nzero_intscore))

  return(results)

  }
