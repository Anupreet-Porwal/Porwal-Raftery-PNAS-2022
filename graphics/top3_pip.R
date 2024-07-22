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
library(car)
library(data.table)
library(openxlsx)
library(SIS)
library(spikeslab)
library(gss)

registerDoParallel(cores=2)

data_with_int_quad <- function(original.data){
  # Standardizing all the x's except the factor variables
  for (i in 1:(ncol(original.data)-1)){
    if(is.factor(original.data[ ,i])==0){
      original.data[ ,i]=(original.data[ ,i]-mean(original.data[ ,i]))/sd(original.data[ ,i])
    }
  }
  
  # centering the Y. Assuming Y is always the last column of original.data
  y.val.0 <- mean(original.data[ , ncol(original.data)])
  original.data[ , ncol(original.data)]=original.data[ , ncol(original.data)]-mean(original.data[ , ncol(original.data)])
  p <- ncol(original.data)-1
  colnames(original.data)[ncol(original.data)]="Y"
  colnames(original.data)[1:p] <- paste("V",1:p,sep="")
  
  
  X.original <- original.data[ ,-ncol(original.data)]
  Y.original <- original.data[ , ncol(original.data)]
  
  colnum <- ncol(X.original)
  coltype <- sapply(X.original, class)
  
  
  # Include all the two way interaction terms; remove the intercept term in next line
  X.original <- model.matrix(Y~.*.,data=original.data)
  X.original <- X.original[ ,-1]
  # for (j in 1:(colnum-1)){
  #   for (k in (j+1):colnum){
  #     X.original <- cbind.data.frame(X.original, X.original[ ,j]*X.original[ ,k])
  #     colnames(X.original)[ncol(X.original)] <-
  #       paste(colnames(X.original)[j],colnames(X.original)[k],sep=".")
  #   }
  #
  # }
  #
  for (j in 1:colnum){
    if(coltype[j]=="numeric"){
      X.original <- cbind.data.frame(X.original, X.original[ ,j]^2)
      colnames(X.original)[ncol(X.original)] <-
        paste(colnames(X.original)[j],2,sep=".")
    }else if(coltype[j]=="factor"){
      X.original[ ,j] <- as.factor(X.original[ ,j])
    }
  }
  
  
  prelim.data <- cbind.data.frame(X.original,Y.original)
  colnames(prelim.data)[ncol(prelim.data)] <- "Y"
  
  return(prelim.data)
  
}


#####
#setwd("C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code")
#source("src/est_varselection_function_011321.R")



#args = commandArgs(TRUE)

datanames <- c("College", "BC-Tmax", "BC-Tmin", "BS-daily",
               "BS-hourly","SML","Diabetes", "Superconductivity", "Ozone", "Boston",
               "Nutrimouse","Multidrug","NIR","Liver")


for(data.num in 1:length(datanames)){

dataname <- datanames[data.num]



if(dataname=="College"){
  data("College")
  prelim.data=College
  
  col.num <- ncol(prelim.data)
  prelim.data=prelim.data[ , c(1,3:(col.num),2)]
  # removed enroll and accept due to causual issues
  prelim.data=subset(prelim.data, select=-c(Enroll,Accept))
  prelim.data$Apps=log(prelim.data$Apps)
  prelim.data$F.Undergrad=log(prelim.data$F.Undergrad)
  prelim.data$P.Undergrad=log(prelim.data$P.Undergrad)
  
}else if(dataname=="BC-Tmax"){
  prelim.data <- read.table("data/BiasCorrection/Bias_correction_ucl.csv",header = TRUE,sep = ",")
  
  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , -c(1:2,colnum)]
  prelim.data <- prelim.data[complete.cases(prelim.data),]
  prelim.data[ , ncol(prelim.data)] <- sqrt(prelim.data[ , ncol(prelim.data)] )
  
  
}else if(dataname=="BC-Tmin"){
  prelim.data <- read.table("data/BiasCorrection/Bias_correction_ucl.csv",header = TRUE,sep = ",")
  
  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , -c(1:2,colnum-1)]
  prelim.data <- prelim.data[complete.cases(prelim.data),]
  prelim.data[ , ncol(prelim.data)] <- sqrt(prelim.data[ , ncol(prelim.data)] )
  
}else if(dataname=="BS-daily"){
  prelim.data=read.table("data/bikesharing/day.csv",header = TRUE,sep = ",")
  
  # remove date and record index since these variables are just indexing variables
  # remove casual and registered since they define the outcome variable of interest
  prelim.data=prelim.data[ ,-c(1,2,14,15)]
  
  # Changing normalised versions to actual numbers based on the read me file
  prelim.data$temp <- prelim.data$temp*41
  prelim.data$atemp <- prelim.data$atemp*50
  prelim.data$windspeed <- prelim.data$windspeed*67
  prelim.data$hum <- prelim.data$hum*100
  
  prelim.data$yr <- factor(format(prelim.data$yr, format="%A"),
                           levels = c("0", "1") , labels = c("2011","2012"))
  prelim.data$weathersit <- factor(format(prelim.data$weathersit, format="%A"),
                                   levels = c("1", "2","3") ,
                                   labels = c("Good","Moderate","Bad"))
  prelim.data$holiday <- factor(format(prelim.data$holiday, format="%A"),
                                levels = c("0", "1") , labels = c("NotHoliDay","Holiday"))
  prelim.data$season <- factor(format(prelim.data$season, format="%A"),
                               levels = c("1", "2","3","4") , labels = c("Spring","Summer","Fall","Winter"))
  prelim.data$mnth <- factor(prelim.data$mnth)
  prelim.data$mnth <- relevel(prelim.data$mnth,ref=6)
  prelim.data$weekday <- factor(prelim.data$weekday)
  prelim.data$weekday <- relevel(prelim.data$weekday,ref=3)
  prelim.data=subset(prelim.data,select = -c(workingday))
  prelim.data$cnt=(prelim.data$cnt)^(1/2)
  
}else if(dataname=="BS-hourly"){
  
  prelim.data=read.table("data/bikesharing/hour.csv",header = TRUE,sep = ",")
  
  # remove date and record index since these variables are just indexing variables
  # remove casual and registered since they define the outcome variable of interest
  prelim.data=prelim.data[ ,-c(1,2,15,16)]
  
  # Changing normalised versions to actual numbers based on the read me file
  prelim.data$temp <- prelim.data$temp*41
  prelim.data$atemp <- prelim.data$atemp*50
  prelim.data$windspeed <- prelim.data$windspeed*67
  prelim.data$hum <- prelim.data$hum*100
  
  prelim.data$yr <- factor(format(prelim.data$yr, format="%A"),
                           levels = c("0", "1") , labels = c("2011","2012"))
  prelim.data$weathersit <- factor(format(prelim.data$weathersit, format="%A"),
                                   levels = c("1", "2","3","4") ,
                                   labels = c("Good","Moderate","Bad","Bad"))
  prelim.data$holiday <- factor(format(prelim.data$holiday, format="%A"),
                                levels = c("0", "1") , labels = c("NotHoliDay","Holiday"))
  prelim.data$season <- factor(format(prelim.data$season, format="%A"),
                               levels = c("1", "2","3","4") , labels = c("Spring","Summer","Fall","Winter"))
  prelim.data$hr[prelim.data$hr %in% 0:5] <- "LateNight"
  prelim.data$hr[prelim.data$hr %in% 6:8] <- "EarlyMorning"
  prelim.data$hr[prelim.data$hr %in% 9:15] <- "Morning"
  prelim.data$hr[prelim.data$hr %in% 16:18] <- "Evening"
  prelim.data$hr[prelim.data$hr %in% 19:23] <- "Night"
  
  prelim.data$hr <- as.factor(prelim.data$hr)
  
  prelim.data$hr <- relevel(prelim.data$hr,ref="Morning")
  
  # Cube root tranformation of dependent variable, season: ref=summer
  # month: ref=may, weekday:ref=wednesday, hr:ref=Morning
  
  prelim.data$season <- relevel(prelim.data$season,ref=2)
  prelim.data$mnth <- factor(prelim.data$mnth)
  prelim.data$mnth <- relevel(prelim.data$mnth,ref=5)
  prelim.data$weekday <- factor(prelim.data$weekday)
  prelim.data$weekday <- relevel(prelim.data$weekday,ref=3)
  prelim.data=subset(prelim.data,select = -c(workingday))
  prelim.data$cnt=(prelim.data$cnt)^(1/3)
  
}else if(dataname=="SML"){
  
  prelim.data <- read.table("data/SML2010/NEW-DATA-2.T15.txt",header = TRUE,sep = " ")
  
  prelim.data <- prelim.data[ , -c(1:2,19:21)]
  colnum <- ncol(prelim.data)
  prelim.data <- prelim.data[ , c(3:colnum,1)]
  prelim.data$X24.Day_Of_Week <- round(prelim.data$X24.Day_Of_Week,digits = 0)
  prelim.data$X24.Day_Of_Week <- factor(prelim.data$X24.Day_Of_Week)
  
}else if(dataname=="Diabetes"){
  data("diabetesI")
  prelim.data <- diabetesI
  
  prelim.data <- cbind.data.frame(prelim.data[ ,-1],prelim.data[ ,1])
  
}else if(dataname=="Superconductivity"){
  prelim.data <- read.table("data/Superconductivity/train.csv",header = TRUE,sep = ",")
  # transforming Y variable with cube root transformation
  prelim.data$critical_temp <- prelim.data$critical_temp^(1/3)
  
}else if(dataname=="Ozone"){
  # Data version used in Liang et Al (2008), Miller 2001
  data("ozone",package = "gss")
  ozone <- ozone[ ,-ncol(ozone)]
  original.data <- cbind.data.frame(ozone[ ,-1], ozone[ ,1])
  
  colnames(original.data)[ncol(original.data)] <- "Ozone"
  original.data$Ozone <- log(original.data$Ozone)
  
  prelim.data <- data_with_int_quad(original.data)
  
  
  
}else if(dataname=="Boston"){
  data("BostonHousing")
  original.data <- BostonHousing
  prelim.data <- data_with_int_quad(original.data)
  
}else if(dataname=="Nutrimouse"){
  library(mixOmics)
  
  data("nutrimouse")
  
  prelim.data <- cbind(nutrimouse$gene,nutrimouse$lipid$C16.0)
  prelim.data <- as.data.frame(prelim.data)
  colnames(prelim.data)[ncol(prelim.data)] <- colnames(nutrimouse$lipid)[2]
  
}else if(dataname=="Multidrug"){
  library(mixOmics)
  data("multidrug")
  
  x.dirty <- multidrug$compound
  x.clean <- x.dirty[ , colSums(is.na(x.dirty))==0]
  
  y <- multidrug$ABC.trans[ ,3]
  
  prelim.data <- cbind(x.clean,y)
  prelim.data <- as.data.frame(prelim.data)
  
}else if(dataname=="NIR"){
  library(chemometrics)
  
  data(NIR)
  # response is square root of glucose
  prelim.data <- cbind(NIR$xNIR,sqrt(NIR$yGlcEtOH$Glucose))
  prelim.data <- as.data.frame(prelim.data)
  
}else if(dataname=="Liver"){
  library(mixOmics)
  data("liver.toxicity")
  
  prelim.data <- cbind(liver.toxicity$gene, liver.toxicity$clinic$Cholesterol.mg.dL.)
  prelim.data <- as.data.frame(prelim.data)
  
}


# Standardizing all the x's except the factor variables
for (i in 1:(ncol(prelim.data)-1)){
  if(is.factor(prelim.data[ ,i])==0){
    prelim.data[ ,i]=(prelim.data[ ,i]-mean(prelim.data[ ,i]))/sd(prelim.data[ ,i])
  }
}

# centering the Y. Assuming Y is always the last column of prelim.data
y.val.0 <- mean(prelim.data[ , ncol(prelim.data)])
prelim.data[ , ncol(prelim.data)]=prelim.data[ , ncol(prelim.data)]-mean(prelim.data[ , ncol(prelim.data)])
p <- ncol(prelim.data)-1
colnames(prelim.data)[1:p] <- paste("V",1:p,sep="")
colnames(prelim.data)[ncol(prelim.data)]="Y"

#p <- ncol(prelim.data)-1
N <- nrow(prelim.data)
#### NO EDITS BELOW THIS LINE (except RESULTS) ####
datamat <- prelim.data

Xmat=model.matrix(Y~., data = datamat)
p <- ncol(Xmat)-1
Y=datamat[ ,ncol(datamat)]

# For high dimension, add the truncated model prior
EBlocal.mod <- bas.lm(Y~.,datamat, prior = "EB-local",method = "MCMC", MCMC.iterations = 10000)#, modelprior = tr.beta.binomial(1,1,N-2)
gsqrtn.mod <- bas.lm(Y~., datamat, prior = "g-prior", alpha = sqrt(N),method = "MCMC", MCMC.iterations = 10000)
hypergprior.mod <- bas.lm(Y~.,datamat, prior = "hyper-g",method = "MCMC", MCMC.iterations = 10000)

coef.mat <- matrix(NA, nrow=p+1, ncol = 3)
colnames(coef.mat) <- c("gsqrtn", "hyperg","EBlocal")

coef.mat[ ,1] <- coef(gsqrtn.mod)$postmean
coef.mat[ ,2] <- coef(hypergprior.mod)$postmean
coef.mat[ ,3] <- coef(EBlocal.mod)$postmean

# pairs(coef.mat, pch = 19, lower.panel = NULL,panel= function(x,y,...){ 
#   points(x,y); abline(a=0, b=1)},main=paste("Pairwise scatter plot for parameter estimates- ",datanames[data.num]))

pip.mat <- matrix(NA, nrow=p+1, ncol = 3)
colnames(pip.mat) <- c("gsqrtn", "hyperg","EBlocal")

pip.mat[ ,1] <- gsqrtn.mod$probne0
pip.mat[ ,2] <- hypergprior.mod$probne0
pip.mat[ ,3] <- EBlocal.mod$probne0

# pairs(pip.mat, pch = 19, lower.panel = NULL,panel= function(x,y,...){ 
#   points(x,y); abline(a=0, b=1)},main=paste("PIPs for",datanames[data.num]))


coef.mat<- data.frame(coef.mat)
pip.mat <- data.frame(pip.mat)


p1 <- ggplot(coef.mat, aes(x=gsqrtn, y=hyperg)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+ 
  xlab("g=sqrt(n)") + ylab("hyper-g")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


p2 <- ggplot(coef.mat, aes(x=gsqrtn, y=EBlocal)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+
  xlab("g=sqrt(n)") + ylab("EB local")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


p3 <- ggplot(coef.mat, aes(x=EBlocal, y=hyperg)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+
  xlab("EB local") + ylab("hyper-g")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))

p4 <- ggplot(pip.mat, aes(x=gsqrtn, y=hyperg)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+
  xlab("g=sqrt(n)") + ylab("hyper-g")+ 
  ylim(0, 1)+
  xlim(0, 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


p5 <- ggplot(pip.mat, aes(x=gsqrtn, y=EBlocal)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+
  xlab("g=sqrt(n)") + ylab("EB local")+ 
  ylim(0, 1)+
  xlim(0, 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


p6 <- ggplot(pip.mat, aes(x=EBlocal, y=hyperg)) + 
  geom_point(shape=21,fill="blue")+
  #geom_smooth(method=lm)+
  geom_abline(slope=1,intercept = 0)+
  xlab("EB local") + ylab("hyper-g") + 
  ylim(0, 1)+
  xlim(0, 1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


library(cowplot)
p.coef <- plot_grid(p1,p2,p3, nrow=1,ncol=3)+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
p.pip <- plot_grid(p4,p5,p6, nrow=1,ncol=3)+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))

p <- plot_grid(p.coef,p.pip, nrow=2,
               labels=c(paste('Estimated coefficients for', 
                              dataname, 'dataset' ),
                        paste('Posterior Inclusion probabilities (PIPs) for', 
                              dataname, 'dataset')), label_size = 8, vjust = 0.1)+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

output_directory <- 'C:/Users/Anupreet Porwal/Dropbox/Research/BMA LASSO/code/BMA code/results/top3_coef_pip/'

pdf( paste(output_directory,dataname,".pdf",sep=""), height=6, width=10)
print(p)
dev.off()

}