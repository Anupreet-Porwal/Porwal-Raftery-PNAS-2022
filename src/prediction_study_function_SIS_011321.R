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


post.samples.gprior<-function(n.post, gprior.obj,gprior.pred.obj, test.set ){
  post.samp<-matrix(NA, nrow = nrow(test.set), ncol = n.post)

  for (w in 1:ncol(post.samp)){
    index<-sample(1:gprior.obj$n.models, 1, replace = TRUE, prob = gprior.obj$postprobs)
    post.mean<-gprior.pred.obj$Ypred[which(gprior.pred.obj$best==index), ]
    mod.sd<-sqrt(gprior.obj$mse[index])
    post.samp[ ,w]<-rnorm(nrow(post.samp), mean=post.mean, sd= mod.sd)
  }
  return(post.samp)
}

residualVariance <- function(coef, x, y)
  sum((y - x %*% coef)^2)/(length(y) - length(coef))

prediction_study=function(Xmat,Y, bootn=100, split=0.75){

  #### Loop over Bootstrap samples #### ####----
  tot.num <- ncol(model.matrix(Y~.,cbind(Xmat,Y)))

  # methods included in the point/ Interval estimation analysis
  methods_Point=c("BMA-bicreg","Spikeslab","UIP", "BIC ","AIC","EB-local","EB-global",
                  "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","SS Lasso", "EMVS", "SCAD",
                  "MCP","LASSO-lambda.min","LASSO-lambda.1se","Elastic Net","JZS","ZS-null","full model")# ,"ScanBMA")
  methods_Interval=c("BMA-bicreg","Spikeslab","UIP", "BIC ","AIC","EB-local","EB-global",
                     "g-sqrtn","g-1","Hyper g","NLP","Horseshoe","JZS","ZS-null","full model")# ,"ScanBMA")

  mse_mat=rmse_mat=mae_mat= r2test.mat=
    matrix(NA,nrow = bootn,ncol = length(methods_Point))
  colnames(mse_mat)=colnames(mae_mat)=colnames(rmse_mat)= colnames(r2test.mat) =methods_Point

  width_mat=coverage_mat=IntervalScore_mat=CRPS_mat=
    matrix(NA, nrow = bootn, ncol = length(methods_Interval)) # EMVS, SS LASSO, LASSO,SCAD,MCP excluded
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

    if(ncol(X.samp)>30){
      comp.time.sis <- system.time({sis.mod <- SIS(X.samp[ ,-1],Y.samp,family = "gaussian")})[3]
      X.samp.scr <- cbind(rep(1,N.train), X.samp[ ,sis.mod$ix+1])
      #X.samp.scr <- as.matrix(X.samp.scr)
      colnames(X.samp.scr) <- c("(Intercept)",colnames(X.samp)[sis.mod$ix+1])
      test.samp.scr <- test.samp[ , c(1, sis.mod$ix+1)]
    }
    #sampdata.scr <- as.data.frame(sampdata.scr)


      #Bicreg BMA
      bma.pred <- matrix(0,nrow = nrow(test), ncol=4)
      if(ncol(X.samp.scr)>30){
        ibma.mod <- iBMA.bicreg(X.samp.scr[ ,-1], Y.samp,OR=20)
        bma.mod <- ibma.mod$bma
        test.samp.scr <- test.samp.scr[ , c("(Intercept)",
                                            colnames(ibma.mod$sortedX)[ibma.mod$currentSet])]
        X.samp.scr <- X.samp.scr[ , c("(Intercept)",colnames(ibma.mod$sortedX)[ibma.mod$currentSet])]
      }else{
        bma.mod <- bicreg(X.samp.scr[ ,-1], Y.samp,OR=100)
      }
      #bma.pred <- predict(bma.mod,test.samp.scr[ ,-1 ], quantiles=c(0.025,0.975))
      colnames(bma.pred) <- c("Predicted","sd","2.5%","97.5%")
      linpred <-  bma.mod$mle %*% t(test.samp.scr)
      predmean <- apply(bma.mod$postprob * linpred, 2, sum)
      bma.mod[["residvar"]] <-  apply(bma.mod$mle, 1, residualVariance, x = X.samp.scr, y = Y.samp)
      objVARterm <- sum(bma.mod$postprob * bma.mod$residvar)
      predSD <- sqrt(objVARterm +
                       apply((predmean - t(linpred))^2, 1,weighted.mean,w=bma.mod$postprob))
      bma.pred[ ,1] <- predmean
      bma.pred[ ,2] <- predSD
      alpha <- 0.025
      for (t in 1:nrow(test.samp)){
        bma.pred[ t,3:4] <-  sapply(c(0.025,0.975), quantBMAnormal,  WEIGHTS=bma.mod$postprob,
                                    MEAN=linpred[,t], SD=predSD[t])
      }
      bma.pred.mat <- post.samples(2500, bma.mod, test.samp.scr)
      bma.phat <- sum(bma.mod$size*bma.mod$postprob)+1 # Include Intercept
      # mse_mat[b,1] <- MSE(y.true.test, bma.pred[ ,1])
      # mae_mat[b,1] <- MAE(y.true.test, bma.pred[ ,1])
      # rmse_mat[b,1] <- rmse(y.true.test, bma.pred[ ,1])
      # r2test.mat[b,1] <- r2test(y.true.test,bma.pred[ ,1])
      #
      # width_mat[b,1] <- mean(bma.pred[ ,4]-bma.pred[ ,3])
      # coverage_mat[b,1] <- sum(((y.true.test>bma.pred[ ,3])&(y.true.test<bma.pred[ ,4])))/length(y.true.test)
      # IntervalScore_mat[b,1] <- meanintscore(bma.pred[ ,4],bma.pred[ ,3], y.true.test, a)
      # CRPS_mat[b,1] <- CRPS_prediction(100, bma.pred.mat ,y.true.test)



      # Spike and slab
      ss.mod <- quiet(lm.spike(Y~.,niter = 4000, data = train, error.distribution = "gaussian"))
      ss.pred <- predict(ss.mod, test[ ,-ncol(test)], burn=500)
      ss.phat <- mean(rowSums(abs(ss.mod$beta)>0))
      ss.conf <- t(apply(ss.pred, 1,quantile, probs=c(0.025,0.975)))


        # Unit information prior
        gpriorUnit.mod <- bas.lm(Y~., train, prior = "g-prior", alpha = N.train,method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        gpriorUnit.pred <- predict(gpriorUnit.mod, test[ ,-ncol(test)], se.fit=TRUE)
        gpriorUnit.pred.mat <- post.samples.gprior(2500, gpriorUnit.mod,gpriorUnit.pred, test.samp)
        gpriorUnit.phat <- sum(gpriorUnit.mod$postprobs*gpriorUnit.mod$size)
        gpriorUnit.conf <- confint(gpriorUnit.pred)

      #BIC g prior
        gpriorBIC.mod <- bas.lm(Y~., train, prior = "BIC",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        gpriorBIC.pred <- predict(gpriorBIC.mod, test[ ,-ncol(test)], se.fit=TRUE)
        gpriorBIC.pred.mat <- post.samples.gprior(2500, gpriorBIC.mod, gpriorBIC.pred, test.samp)
        gpriorBIC.phat <- sum(gpriorBIC.mod$postprobs*gpriorBIC.mod$size)
        gpriorbic.conf <- confint(gpriorBIC.pred)


      #AIC
        AIC.mod <- bas.lm(Y~., train, prior = "AIC",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        AIC.pred <- predict(AIC.mod, test[ ,-ncol(test)], se.fit=TRUE)
        AIC.pred.mat <- post.samples.gprior(2500, AIC.mod,AIC.pred, test.samp)
        AIC.phat <- sum(AIC.mod$postprobs*AIC.mod$size)
        AIC.conf <- confint(AIC.pred)


      #EB-local
        EBlocal.mod <- bas.lm(Y~., train, prior = "EB-local",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        EBlocal.pred <- predict(EBlocal.mod, test[ ,-ncol(test)], se.fit=TRUE)
        EBlocal.pred.mat <- post.samples.gprior(2500, EBlocal.mod,EBlocal.pred, test.samp)
        EBlocal.phat <- sum(EBlocal.mod$postprobs*EBlocal.mod$size)
        EBlocal.conf <- confint(EBlocal.pred)


      #EB-global
        EBglobal.mod <- bas.lm(Y~., train, prior = "EB-global",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        EBglobal.pred <- predict(EBglobal.mod, test[ ,-ncol(test)], se.fit=TRUE)
        EBglobal.pred.mat <- post.samples.gprior(2500, EBglobal.mod,EBglobal.pred, test.samp)
        EBglobal.phat <- sum(EBglobal.mod$postprobs*EBglobal.mod$size)
        EBglobal.conf <- confint(EBglobal.pred)


      #g = sqrt(N)
        gsqrtn.mod <- bas.lm(Y~., train, prior = "g-prior", alpha = sqrt(N.train),method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        gsqrtn.pred <- predict(gsqrtn.mod, test[ ,-ncol(test)], se.fit=TRUE)
        gsqrtn.pred.mat <- post.samples.gprior(2500, gsqrtn.mod,gsqrtn.pred, test.samp)
        gsqrtn.phat <- sum(gsqrtn.mod$postprobs*gsqrtn.mod$size)
        gsqrtn.conf <- confint(gsqrtn.pred)


      #g = 1
        g1.mod <- bas.lm(Y~., train, prior = "g-prior", alpha = 1,method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        g1.pred <- predict(g1.mod, test[ ,-ncol(test)], se.fit=TRUE)
        g1.pred.mat <- post.samples.gprior(2500, g1.mod,g1.pred, test.samp)
        g1.phat <- sum(g1.mod$postprobs*g1.mod$size)
        g1.conf <- confint(g1.pred)


      # Hyper g priors
        hypergprior.mod <- bas.lm(Y~.,train, prior = "hyper-g",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
        hypergprior.pred <- predict(hypergprior.mod, test[ ,-ncol(test)], se.fit=TRUE)
        hypergprior.pred.mat <- post.samples.gprior(2500, hypergprior.mod, hypergprior.pred, test.samp)
        hyperg.phat <- sum(hypergprior.mod$postprobs*hypergprior.mod$size)
        hypergprior.conf <- confint(hypergprior.pred)

        #Non-Local Priors
        nlp.pred <- matrix(0,nrow = nrow(test), ncol=3)
          nlp.modselect <- quiet(modelSelection(y=Y.samp,x=X.samp, data= train, scale = FALSE,center = FALSE))
          nlp.postsamp <- rnlp(msfit=nlp.modselect, niter = 2500,burnin = 1000)
          nlp.betasamp <- nlp.postsamp[ , 2:(ncol(nlp.postsamp)-1)]
          nlp.phisamp <- nlp.postsamp[ , ncol(nlp.postsamp)]
          nlp.pred.mean <- test.samp %*% t(nlp.betasamp)
          nlp.pred.mat <- t(apply(nlp.pred.mean, 1, rnorm, n=ncol(nlp.pred.mean), sd=sqrt(nlp.phisamp)))
          colnames(nlp.pred) <- c("Predicted", "2.5%", "97.5%")
          nlp.pred[, 2:3] <-  t(apply(nlp.pred.mat, 1, quantile, probs=c(0.025,0.975)))
          nlp.pred[, 1] <-  rowMeans(nlp.pred.mat)
          nlp.phat <- mean(rowSums(nlp.modselect$postSample))



        #Horeshoe
        HS.pred <- matrix(0,nrow = nrow(test), ncol=3)
        HS.mod <- quiet(horseshoe(Y.samp,X.samp,method.tau = "halfCauchy", method.sigma = "Jeffreys",nmc = 2500))
        HS.vs <- HS.var.select(HS.mod, Y.samp, method = "intervals")
        colnames(HS.pred) <- c("Predicted", "2.5%", "97.5%")
        HS.pred.mat.mean <-  test.samp %*% HS.mod$BetaSamples
        HS.pred.mat <- t(apply(HS.pred.mat.mean, 1, rnorm, n=ncol(HS.pred.mat.mean), sd=sqrt(HS.mod$Sigma2Samples)))
        HS.pred[, 2:3] <-  t(apply(HS.pred.mat, 1, quantile, probs=c(0.025,0.975)))
        HS.pred[, 1] <-  test.samp %*% HS.mod$BetaHat
          alpha <- 0.5
          effsamp <- ncol(HS.mod$BetaSamples)
          left.50 <- floor(alpha*effsamp/2)
          right.50 <- ceiling((1-alpha/2)*effsamp)
          BetaSort <- apply(HS.mod$BetaSamples, 1, sort, decreasing = F)
          left.points <- BetaSort[left.50, ]
          right.points <- BetaSort[right.50, ]
          HS.vs5 <- as.numeric( 1 - ( (left.points <= 0) & (right.points >= 0) ))
          #Note sure how to calculate ; currently doing 50% credible intervals
          HS.phat <- sum(HS.vs5)


        #Spike & slab LASSO
          sslasso.mod <- SSLASSO(X.samp[ ,-1], Y.samp, variance = "unknown") # you dont need to give a coefficient column seperately
          coef.sslasso <- sslasso.mod$beta[ ,100]
          sslasso.pred <- rep(sslasso.mod$intercept[ ,100], nrow(test.samp))+test.samp[ ,-1] %*% coef.sslasso
          sslasso.phat <- sum(abs(coef.sslasso)>0)+1 # Include Intercept




        #EMVS
          p <- ncol(X.samp)-1
          epsilon <- 10^{-5}
          a.em <- b.em <- 1
          beta.init <- rep(1,p)
          sigma.init <- 1
          v0 <- seq(0.001,15,length.out=1000)
          v1 <- 1000
          # independent =TRUE, temp=1 and betainit set to default
          EMVS.mod <- quiet(EMVS(Y= Y.samp, X = X.samp[ ,-1],v0=v0,v1=v1,
                                 type="betabinomial",
                                 sigma_init = sigma.init ,
                                 epsilon=epsilon,a=a.em,b=b.em))
          # EMVS.mod <- quiet(EMVS(Y= Y.samp, X = X.samp[ ,-1],v0=v0,v1=v1,
          #                        type="betabinomial",beta_init = beta.init,
          #                        sigma_init = sigma.init ,independent = FALSE,
          #                        epsilon=epsilon,a=a.em,b=b.em, temperature=0.2))
          #EMVS.best <- quiet(EMVSbest(EMVS.mod))
          emvs.coef <- EMVS.mod$betas[1, ]
          emvs.pred <- test.samp[ ,-1] %*% emvs.coef
          EMVS.phat <- length(which(EMVS.mod$prob_inclusion[1,]>=0.5))

          #SCAD
            scad.mod <-  cv.ncvreg(X.samp[ ,-1], Y.samp, family = "gaussian", penalty = "SCAD")
            scad.pred <- predict(scad.mod, test.samp[ ,-1])
            scad.phat <- sum(abs(coef(scad.mod))>0)

          #MCP
            mcp.mod <-  cv.ncvreg(X.samp[ ,-1], Y.samp, family = "gaussian", penalty = "MCP")
            mcp.pred <- predict(mcp.mod, test.samp[ ,-1])
            mcp.phat <- sum(abs(coef(mcp.mod))>0)

          # LASSO lambda min
            lasso.mod <- cv.glmnet(X.samp[ ,-1], Y.samp, family = "gaussian")
            lasso.pred.min <- predict(lasso.mod, test.samp[ ,-1], s="lambda.min")
            lasso.min.phat <- sum(abs(coef(lasso.mod, s="lambda.min"))>0)

            # Lasso Lambda 1se
            lasso.pred.1se <- predict(lasso.mod, test.samp[ ,-1], s="lambda.1se")
            lasso.1se.phat <- sum(abs(coef(lasso.mod, s="lambda.1se"))>0)


          # Elastic net
            alp <- seq(0.1, 0.95, 0.05)
            search <- foreach(i = alp, .combine = rbind, .packages = "glmnet") %dopar% {
              cv <- cv.glmnet(X.samp[ , -1], Y.samp, family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
              data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
            }
            cv3 <- search[search$cvm == min(search$cvm), ]
            elastic.mod <- glmnet(X.samp[ , -1], Y.samp, family = "gaussian", lambda = cv3$lambda.min, alpha = cv3$alpha)
            elastic.pred <- predict(elastic.mod, test.samp[ ,-1])
            elastic.phat <- sum(abs(coef(elastic.mod))>0)

        #JZS
      JZS.mod <- bas.lm(Y~., train, prior = "JZS",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
      JZS.pred <- predict(JZS.mod, test[ ,-ncol(test)], se.fit=TRUE)
      JZS.pred.mat <- post.samples.gprior(2500, JZS.mod,JZS.pred, test.samp)
      JZS.phat <- sum(JZS.mod$postprobs*JZS.mod$size)
      JZS.conf <- confint(JZS.pred)



    #ZS-null
      ZSnull.mod <- bas.lm(Y~., train, prior = "ZS-null",method = "MCMC",MCMC.iterations = 10000, modelprior = tr.beta.binomial(1,1,N.train-2))
      ZSnull.pred <- predict(ZSnull.mod, test[ ,-ncol(test)], se.fit=TRUE)
      ZSnull.pred.mat <- post.samples.gprior(2500, ZSnull.mod,ZSnull.pred, test.samp)
      ZSnull.phat <- sum(ZSnull.mod$postprobs*ZSnull.mod$size)
      ZSnull.conf <- confint(ZSnull.pred)



    full.mod.pred <- matrix(NA, nrow = N.test, ncol = 3)
    fullmod.conf <- full.mod.pred[ ,2:3]
    full.phat <- NA
    # full model
    if(p<N.train){
        full.mod <- lm(Y~.,data = train)
        full.mod.pred <- predict(full.mod, newdata=test[ ,-ncol(test)], interval="prediction")
        full.phat <- length(full.mod$coefficients)
        fullmod.conf <- full.mod.pred[ ,2:3]
    }

    pred.mat <- list(bma.pred[ ,1],rowMeans(ss.pred), gpriorUnit.pred$fit,gpriorBIC.pred$fit,
                     AIC.pred$fit, EBlocal.pred$fit, EBglobal.pred$fit, gsqrtn.pred$fit,
                     g1.pred$fit, hypergprior.pred$fit, nlp.pred[ ,1], HS.pred[ ,1],
                     sslasso.pred, emvs.pred[ ,1], scad.pred, mcp.pred,lasso.pred.min,
                     lasso.pred.1se,elastic.pred ,JZS.pred$fit, ZSnull.pred$fit, full.mod.pred[ ,1])
    names(pred.mat) <- methods_Point

    for (i in 1:length(methods_Point)){
      mse_mat[b,i] <- MSE(y.true.test, pred.mat[[i]])
      rmse_mat[b,i] <- rmse(y.true.test, pred.mat[[i]])
      mae_mat[b,i] <- MAE(y.true.test, pred.mat[[i]])
      r2test.mat[b,i] <- r2test(y.true.test,pred.mat[[i]])
    }

    phat.mat[b, ] <- cbind(bma.phat, ss.phat, gpriorUnit.phat,gpriorBIC.phat, AIC.phat,
                           EBlocal.phat,EBglobal.phat,gsqrtn.phat, g1.phat, hyperg.phat,
                           nlp.phat,HS.phat,sslasso.phat, EMVS.phat, scad.phat, mcp.phat,
                           lasso.min.phat,lasso.1se.phat, elastic.phat, JZS.phat, ZSnull.phat,
                           full.phat)

    conf.mat <- list(bma.pred[ ,3:4], ss.conf, gpriorUnit.conf[ ,1:2], gpriorbic.conf[ ,1:2],
                     AIC.conf[ ,1:2],EBlocal.conf[ ,1:2],EBglobal.conf[ ,1:2], gsqrtn.conf[ ,1:2],
                     g1.conf, hypergprior.conf[ ,1:2],nlp.pred[,2:3], HS.pred[ ,2:3],
                     JZS.conf[ ,1:2],ZSnull.conf[ ,1:2],fullmod.conf)

    for(j in 1:length(methods_Interval)){
      conf.method <- conf.mat[[j]]
      width_mat[b,j] <- mean(conf.method[ ,2]-conf.method[ ,1])
      coverage_mat[b,j] <- sum(((y.true.test>conf.method[ ,1])&(y.true.test<conf.method[ ,2])))/ length(y.true.test)
      IntervalScore_mat[b,j] <- meanintscore(conf.method[ ,2],conf.method[ ,1], y.true.test, a)
    }

    methods.pred.mat <- list(bma.pred.mat ,ss.pred, gpriorBIC.pred.mat, gpriorBIC.pred.mat,
                             AIC.pred.mat,EBlocal.pred.mat,EBglobal.pred.mat,gsqrtn.pred.mat,
                             g1.pred.mat,hypergprior.pred.mat,nlp.pred.mat,HS.pred.mat,
                             JZS.pred.mat,ZSnull.pred.mat)

    for(j in 1:(length(methods_Interval)-1)){
      CRPS_mat[b,j] <- CRPS_prediction(100, methods.pred.mat[[j]]  ,y.true.test)

    }

  }



  #### Summary plots and tables #### ####----
  # par(mfrow=c(2,4))
  #
  # for (w in 1:ncol(width_mat)){
  #   hist(width_mat[ ,w], main = paste("PI-",colnames(width_mat)[w]), xlab = "Avg Width")
  # }


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
