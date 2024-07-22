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

pos_vector <- rep(3, nrow(dat.dim))
pos_vector[rownames(dat.dim) %in% c("BC-Tmin", "College","NIR","Superconductivity")] <- 1

plot(n~p, data = dat.dim, log="xy",pch=19, cex=1,cex.lab=1.5,cex.axis=1.5,col="blue")
abline(0,1, col="red", lwd=3)
text(n~p, labels=rownames(dat.dim),data=dat.dim, cex=0.7, font=0,pos=pos_vector)
 

###
plot((100-p0*100/p)~log(p), col="blue", pch=19, cex=1,data=dat.dim, ylab="Zero coefficients(%)",cex.lab=1.5,cex.axis=1.5)

abline(h=80,col="red", lty=2)
abline(h=40,col="red",lty=2)


plot(R2~log(p), col="blue", pch=19, cex=1,data=dat.dim)











