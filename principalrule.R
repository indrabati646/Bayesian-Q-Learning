#Principal treatment regime for Simulation Settings 1, 2 and 3

source("/Users/indrabatibhattacharya/Documents/Compliance/metropolisNoTrans.R")
source("/Users/indrabatibhattacharya/Documents/QL_code/qsim2stage.R")
source("/Users/indrabatibhattacharya/Documents/Compliance/metropolisC.R")
source("/Users/indrabatibhattacharya/Documents/QL_code/metropolisy.R")
library(dirichletprocess)
library(mnormt)
library(rootSolve)
library(matrixcalc)
library(condMVNorm)
library(KScorrect)
library(MASS)
library(mvtnorm)
library(ggplot2)


#Grid for the compliance values for the heatmaps of the posterior optimality probability
D1su <- D2su <- seq(0,1,0.05); D21s <- D22s <- seq(0,1,0.05)
Ds <- expand.grid(D1su=D1su,D2su=D2su)

#Function to generate the heatmaps of the posterior optimality probability
qlfun_pri <- function(data){
  data <- matrix(unlist(data),nrow=nrow(data),ncol=ncol(data))
  colnames(data) <- c("a1","a2","D1","D2","D3","D4","s","y","y1","X1","X2","X3","u")
  
  ind01 <- which(data[,"a1"]==1 & data[,"s"]==1)
  ind02 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1)
  ind03 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1)
  ind04 <- which(data[,"a1"]==-1 & data[,"s"]==1)
  ind05 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1)
  ind06 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1)
  
  y1 <- data[which(data[,"a1"]==1 & data[,"s"]==1),]
  y2 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1),]
  y3 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1),]
  y4 <- data[which(data[,"a1"]==-1 & data[,"s"]==1),]
  y5 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1),]
  y6 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1),]
  
  x1 <- data[,"X1"]; x2 <- data[,"X2"]; x3 <- data[,"X3"]
  x <- cbind(x1,x2,x3)
  a1<-data[,"a1"]; a2 <- data[,"a2"]
  
  xx1 <- y1[,c("X1","X2","X3")]
  xx2 <- y2[,c("X1","X2","X3")]
  xx3 <- y3[,c("X1","X2","X3")]
  xx4 <- y4[,c("X1","X2","X3")]
  xx5 <- y5[,c("X1","X2","X3")]
  xx6 <- y6[,c("X1","X2","X3")]
  
  aa1<-y1[,"a1"]; aa2<-y2[,"a1"]; aa3<-y3[,"a1"]
  aa4<-y4[,"a1"]; aa5<-y5[,"a1"]; aa6<-y6[,"a1"]
  aa21<-y1[,"a2"]; aa22<-y2[,"a2"]; aa23<-y3[,"a2"]
  aa24<-y4[,"a2"]; aa25<-y5[,"a2"]; aa26<-y6[,"a2"]
  
  n1 <- nrow(y1); n2 <- nrow(y2); n3 <- nrow(y3); n4 <- nrow(y4) 
  n5 <- nrow(y5); n6 <- nrow(y6)
  n <- n1+n2+n3+n4+n5+n6 
  y <- data[,"y"]
  yy1 <- y1[,"y"]; yy2 <- y2[,"y"]; yy3 <- y3[,"y"]; yy4 <- y4[,"y"]
  yy5 <- y5[,"y"]; yy6 <- y6[,"y"]
  data <- dt <- rbind(y1,y2,y3,y4,y5,y6)
  x <- dt[,c("X1","X2","X3")]
  
  DD1 <- cbind(y1[,"D1"],0,0,0)
  DD2 <- cbind(y2[,"D1"],0,y2[,"D3"],0)
  DD3 <- cbind(y3[,"D1"],0,0,y3[,"D4"])
  DD4 <- cbind(0,y4[,"D2"],0,0)
  DD5 <- cbind(0,y5[,"D2"],y5[,"D3"],0)
  DD6 <- cbind(0,y6[,"D2"],0,y6[,"D4"])
  DD1[,2:4] <- cbind(rtruncnorm(n1,mean(DD1[,1]),0.1,a=0,b=1),0,0)
  DD2[,c(2,4)] <- cbind(rtruncnorm(n2,mean(DD2[,1]),0.1,a=0,b=1),
                        rtruncnorm(n2,mean(DD2[,3]),0.1,a=0,b=1))
  DD3[,2:3] <- cbind(rtruncnorm(n3,mean(DD3[,1]),0.1,a=0,b=1),
                     rtruncnorm(n3,mean(DD3[,4]),0.1,a=0,b=1))
  DD4[,c(1,3,4)] <-cbind(rtruncnorm(n4,mean(DD4[,2]),0.1,a=0,b=1),0,0)
  DD5[,c(1,4)] <- cbind(rtruncnorm(n5,mean(DD5[,2]),0.1,a=0,b=1),
                        rtruncnorm(n5,mean(DD5[,3]),0.1,a=0,b=1))
  DD6[,c(1,3)] <- cbind(rtruncnorm(n6,mean(DD6[,2]),0.1,a=0,b=1),
                        rtruncnorm(n6,mean(DD6[,4]),0.1,a=0,b=1))
  D <- rbind(DD1,DD2,DD3,DD4,DD5,DD6)
  h <- rbind(DD2,DD3,DD5,DD6); h2<- DD2; h3<-DD3; h5<- DD5; h6<- DD6
  
  
  
  K <-6; mcmc <- 20000; P<-2
  dim.cov <- 2;dimx <- dim.cov+1
  R<- matrix(c(1,.2,.2,.2,.2,1,.2,.2,.2,.2,1,.2,.2,.2,.2,1),nrow=4,ncol=4)
  yyi1 <- yym1 <- matrix(0,n2+n3,mcmc-2)
  yyi2 <- yym2 <- matrix(0,n5+n6,mcmc-2)
  
  
  para.dd1 <- matrix(nrow=mcmc,ncol=K+K+2+K+5+1+n2+n3+n5+n6)
  para.dd1[1,] <-para.dd1[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[(n1+1):(n1+n2+n3),1]~rbind(xx2,xx3)[,1]+
                                                                    rbind(xx2,xx3)[,2]))[1],K),
                                   coefficients(lm(D[(n1+1):(n1+n2+n3),1]~rbind(xx2,xx3)[,1]+
                                                     rbind(xx2,xx3)[,2]))[-1],
                                   rep(2,K),1,1,1,mean(D[(n1+1):(n1+n2+n3),1]),1,1, 
                                   sample(1:K,n2+n3+n5+n6, replace=TRUE))
  
  para.dd2 <- matrix(nrow=mcmc,ncol=K+K+2+K+5+1+n2+n3+n5+n6)
  para.dd2[1,] <-para.dd2[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[(n1+n2+n3+n4+1):n,2]~rbind(xx5,xx6)[,1]+
                                                                    rbind(xx5,xx6)[,2]))[1],K),
                                   coefficients(lm(D[(n1+n2+n3+n4+1):n,2]~rbind(xx5,xx6)[,1]+
                                                     rbind(xx5,xx6)[,2]))[-1],
                                   rep(2,K),1,1,1,mean(D[(n1+n2+n3+n4+1):n,2]),1,1, 
                                   sample(1:K,n2+n3+n5+n6, replace=TRUE))
  
  
  para.dd3 <- matrix(nrow=mcmc,ncol=K+K+3+K+5+1+n2+n3+n5+n6)
  para.dd3[1,] <-para.dd3[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[c(ind02,ind05),3]~rbind(xx2,xx5)[,1]+rbind(xx2,xx5)[,2]+
                                                                    rbind(xx2,xx5)[,3]))[1],K),
                                   coefficients(lm(D[c(ind02,ind05),3]~rbind(xx2,xx5)[,1]+rbind(xx2,xx5)[,2]+rbind(xx2,xx5)[,3]))[-1],
                                   rep(2,K),1,1,1,mean(rbind(DD2,DD5)[,3]),1,1, 
                                   sample(1:K, n2+n3+n5+n6, replace=TRUE))
  
  para.dd4 <- matrix(nrow=mcmc,ncol=K+K+3+K+5+1+n2+n3+n5+n6)
  para.dd4[1,] <-para.dd4[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[c(ind03,ind06),4]~rbind(xx3,xx6)[,1]+
                                                                    rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]))[1],K),
                                   coefficients(lm(D[c(ind03,ind06),4]~rbind(xx3,xx6)[,1]+rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]))[-1],
                                   rep(2,K),1,1,1,mean(rbind(DD3,DD6)[,4]),1,1, 
                                   sample(1:K, n2+n3+n5+n6, replace=TRUE))
  
  para.C<-matrix(nrow = mcmc, ncol = 12)
  para.C[1,] <- para.C[2,] <- c(rep(NA,6),rep(0.1,6))
  
  ind1 <- (K+1):(2*K) # intercept
  ind3<- (2*K+1):(2*K+2)
  ind2 <- (2*K+3):(3*K+2) # variance
  ind4 <- 3*K+3 # alpha
  ind5 <- ind4+1 # alpha_beta0
  ind6 <- ind5+1 # alpha_sigma
  ind7 <- ind6+1 # mu_beta0
  ind8 <- ind7+1 # sigma_beta0
  
  ind9 <- length(para.dd1[1,])
  
  ind13<- (2*K+1):(2*K+3)
  ind12 <- (2*K+4):(3*K+3) # variance
  ind14 <- 3*K+4 # alpha
  ind15 <- ind14+1 # alpha_beta0
  ind16 <- ind15+1 # alpha_sigma
  ind17 <- ind16+1 # mu_beta0
  ind18 <- ind17+1 # sigma_beta0
  
  ind19 <- length(para.dd3[1,])
  
  
  cov.d11 <- diag(vcov(lm(D[(n1+1):(n1+n2+n3),1]~rbind(xx2,xx3)[,1]+
                            rbind(xx2,xx3)[,2]))[1,1],K)
  cov.d12 <- diag(vcov(lm(D[(n1+n2+n3+n4+1):n,2]~rbind(xx5,xx6)[,1]+
                            rbind(xx5,xx6)[,2]))[1,1],K)
  cov.d13 <- diag(vcov(lm(rbind(DD2,DD5)[,3]~rbind(xx2,xx5)[,1]+rbind(xx2,xx5)[,2]+
                            rbind(xx2,xx5)[,3]))[1,1],K)
  cov.d14 <- diag(vcov(lm(rbind(DD3,DD6)[,4]~rbind(xx3,xx6)[,1]+rbind(xx3,xx6)[,2]+
                            rbind(xx3,xx6)[,3]))[1,1],K)
  
  cov2.d11 <- diag(diag(vcov(lm(D[(n1+1):(n1+n2+n3),1]~rbind(xx2,xx3)[,1]+
                                  rbind(xx2,xx3)[,2]))[-1,-1]),2)
  cov2.d12 <- diag(diag(vcov(lm(D[(n1+n2+n3+n4+1):n,2]~rbind(xx5,xx6)[,1]+
                                  rbind(xx5,xx6)[,2]))[-1,-1]),2)
  cov2.d13 <- diag(diag(vcov(lm(rbind(DD2,DD5)[,3]~rbind(xx2,xx5)[,1]+rbind(xx2,xx5)[,2]+
                                  rbind(xx2,xx5)[,3]))[-1,-1]),3)
  cov2.d14 <- diag(diag(vcov(lm(rbind(DD3,DD6)[,4]~rbind(xx3,xx6)[,1]+rbind(xx3,xx6)[,2]+
                                  rbind(xx3,xx6)[,3]))[-1,-1]),3)
  
  
  K1 <- 6
  para.y1 <- matrix(nrow=mcmc,ncol=K1+K1+8+K1+5+1+n2+n5)
  para.y1[1,] <-para.y1[2,] <- c(rep(1/K1,K1),rep(coefficients(lm(c(yy2,yy5)~rbind(xx2,xx5)[,1]+
                                                                    rbind(xx2,xx5)[,2]+rbind(xx2,xx5)[,3]+rbind(DD2,DD5)[,1]+
                                                                    rbind(DD2,DD5)[,2]+rbind(DD2,DD5)[,3]+rbind(DD2,DD5)[,4]+
                                                                    c(aa2,aa5)))[1],K1),
                                 coefficients(lm(c(yy2,yy5)~rbind(xx2,xx5)[,1]+
                                                   rbind(xx2,xx5)[,2]+rbind(xx2,xx5)[,3]+rbind(DD2,DD5)[,1]+
                                                   rbind(DD2,DD5)[,2]+rbind(DD2,DD5)[,3]+rbind(DD2,DD5)[,4]+
                                                   c(aa2,aa5)))[-1],
                                 rep(2,K1),1,1,1,mean(c(yy2,yy5)),1,1, 
                                 sample(1:K1,n2+n5, replace=TRUE))
  
  para.y2 <- matrix(nrow=mcmc,ncol=K1+K1+8+K1+5+1+n3+n6)
  para.y2[1,] <-para.y2[2,] <- c(rep(1/K1,K1),rep(coefficients(lm(c(yy3,yy6)~rbind(xx3,xx6)[,1]+
                                                                    rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]+rbind(DD3,DD6)[,1]+
                                                                    rbind(DD3,DD6)[,2]+rbind(DD3,DD6)[,3]+rbind(DD3,DD6)[,4]+c(aa3,aa6)))[1],K1),
                                 coefficients(lm(c(yy3,yy6)~rbind(xx3,xx6)[,1]+
                                                   rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]+rbind(DD3,DD6)[,1]+
                                                   rbind(DD3,DD6)[,2]+rbind(DD3,DD6)[,3]+rbind(DD3,DD6)[,4]+c(aa3,aa6)))[-1],
                                 rep(2,K1),1,1,1,mean(c(yy3,yy6)),1,1, 
                                 sample(1:K1,n3+n6, replace=TRUE))
  
  
  indy1 <- (K1+1):(2*K1) # intercept
  indy3<- (2*K1+1):(2*K1+8)
  indy2 <- (2*K1+9):(3*K1+8) # variance
  indy4 <- 3*K1+9 # alpha
  indy5 <- indy4+1 # alpha_beta0
  indy6 <- indy5+1 # alpha_sigma
  indy7 <- indy6+1 # mu_beta0
  indy8 <- indy7+1 # sigma_beta0
  
  indy91 <- length(para.y1[1,]); indy92 <- length(para.y2[1,])
  
  
  cov.y1 <- diag(vcov(lm(c(yy2,yy5)~rbind(xx2,xx5)[,1]+rbind(xx2,xx5)[,2]+rbind(xx2,xx5)[,3]+
                           rbind(DD2,DD5)[,1]+rbind(DD2,DD5)[,2]+
                           rbind(DD2,DD5)[,3]+rbind(DD2,DD5)[,4]+c(aa2,aa5)))[1,1],K1)
  
  cov2.y1 <- diag(diag(vcov(lm(c(yy2,yy5)~rbind(xx2,xx5)[,1]+
                                 rbind(xx2,xx5)[,2]+rbind(xx2,xx5)[,3]+
                                 rbind(DD2,DD5)[,1]+rbind(DD2,DD5)[,2]+
                                 rbind(DD2,DD5)[,3]+rbind(DD2,DD5)[,4]+c(aa2,aa5)))[-1,-1]),8)
  
  cov.y2 <- diag(vcov(lm(c(yy3,yy6)~rbind(xx3,xx6)[,1]+rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]+
                           rbind(DD3,DD6)[,1]+rbind(DD3,DD6)[,2]+
                           rbind(DD3,DD6)[,3]+rbind(DD3,DD6)[,4]+c(aa3,aa6)))[1,1],K1)
  
  cov2.y2 <- diag(diag(vcov(lm(c(yy3,yy6)~rbind(xx3,xx6)[,1]+
                                 rbind(xx3,xx6)[,2]+rbind(xx3,xx6)[,3]+
                                 rbind(DD3,DD6)[,1]+rbind(DD3,DD6)[,2]+
                                 rbind(DD3,DD6)[,3]+rbind(DD3,DD6)[,4]+c(aa3,aa6)))[-1,-1]),8)
  acc2 <-acc3 <-acc5<-acc6<-0
  cc <-0
  ME <- list()
  psu <- c()
  psu[1:n1] <- y1[,"y"]
  psu1 <- matrix(0,n2+n3,mcmc-2)
  psu[(n1+n2+n3+1):(n1+n2+n3+n4)] <-  y4[,"y"]
  x56 <- rbind(xx5,xx6)
  psu2 <- matrix(0, n5+n6,mcmc-2)
  for(t in 3:mcmc){
    print(t)
    if(t > 200){
      cov.y1 <- cov(para.y1[3:(t-1),(K1+1):(2*K1)])
      cov2.y1 <- cov(para.y1[3:(t-1),indy3])
      cov.y2 <- cov(para.y2[3:(t-1),(K1+1):(2*K1)])
      cov2.y2 <- cov(para.y2[3:(t-1),indy3])
    }
    para.y1[t,] <- metropolisy(Y=c(yy2,yy5), X=rbind(cbind(xx2,DD2,aa2),
                                                     cbind(xx5,DD5,aa5)),
                               w_pre=para.y1[t-1,1:K1],
                               beta0_pre=para.y1[t-1,indy1], 
                               beta_pre=para.y1[t-1,indy3],
                               sigma_pre=para.y1[t-1,indy2],
                               alpha_pre=para.y1[t-1,indy4], 
                               alpha_beta0_pre=para.y1[t-1,indy5], 
                               alpha_sigma_pre=para.y1[t-1,indy6],
                               mu_beta0_pre=para.y1[t-1,indy7], 
                               sigma_beta0_pre=para.y1[t-1,indy8], K=K1, 
                               eps=0.1,del1=15, del2=10, del3=20, 
                               cov1=cov.y1,cov2=cov2.y1,gamma=para.y1[t-1,(indy91-n2-n5)], 
                               Z=para.y1[t-1,(indy91-n2-n5+1):indy91])
    
    para.y2[t,] <- metropolisy(Y=c(yy3,yy6), X=rbind(cbind(xx3,DD3,aa3),cbind(xx6,DD6,aa6)),
                               w_pre=para.y2[t-1,1:K1],
                               beta0_pre=para.y2[t-1,indy1], 
                               beta_pre=para.y2[t-1,indy3],
                               sigma_pre=para.y2[t-1,indy2],
                               alpha_pre=para.y2[t-1,indy4], 
                               alpha_beta0_pre=para.y2[t-1,indy5], 
                               alpha_sigma_pre=para.y2[t-1,indy6],
                               mu_beta0_pre=para.y2[t-1,indy7], 
                               sigma_beta0_pre=para.y2[t-1,indy8], K=K1, 
                               eps=0.1,del1=15, del2=10, del3=20, 
                               cov1=cov.y2,cov2=cov2.y2,gamma=para.y2[t-1,(indy92-n3-n6)], 
                               Z=para.y2[t-1,(indy92-n3-n6+1):indy92])
    
    if(t > 200){
      cov.d11 <- cov(para.dd1[3:(t-1),(K+1):(2*K)])
      cov2.d11 <- cov(para.dd1[3:(t-1),ind3])
    }
    h <- ifelse(is.nan(h), -3, h)
    para.dd1[t,] <- metropolis(h=h, Y=rbind(DD2,DD3,DD5,DD6)[,1], X=rbind(xx2,xx3,xx5,xx6)[,-3],
                               R=R, w_pre=para.dd1[t-1,1:K],
                               beta0_pre=para.dd1[t-1,ind1], 
                               beta_pre=para.dd1[t-1,ind3],
                               sigma_pre=para.dd1[t-1,ind2],
                               alpha_pre=para.dd1[t-1,ind4], 
                               alpha_beta0_pre=para.dd1[t-1,ind5], 
                               alpha_sigma_pre=para.dd1[t-1,ind6],
                               mu_beta0_pre=para.dd1[t-1,ind7], 
                               sigma_beta0_pre=para.dd1[t-1,ind8], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=1, 
                               cov1=cov.d11,cov2=cov2.d11,gamma=para.dd1[t-1,(ind9-n2-n3-n5-n6)], 
                               Z=para.dd1[t-1,(ind9-n2-n3-n5-n6+1):ind9])
    
    Zm1 <- para.dd1[t,(ind9-n2-n3-n5-n6+1):ind9][1:(n2+n3)]
    h[1:(n2+n3),1] <- qnorm(pmin(1-0.1^10,pmax(0.1^10, ptruncnorm(rbind(DD2,DD3)[,1],a=0,b=1,
                                                                  para.dd1[t,ind1][Zm1]+
                                                                    rbind(xx2,xx3)[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                                                  sqrt(para.dd1[t,ind2])[Zm1]))))
    
    
    
    if(t > 200){
      cov.d12 <- cov(para.dd2[3:(t-1),(K+1):(2*K)])
      cov2.d12 <- cov(para.dd2[3:(t-1),ind3])
    }
    
    h <- ifelse(is.nan(h), -3, h)
    para.dd2[t,] <- metropolis(h=h, Y=rbind(DD2,DD3,DD5,DD6)[,2], 
                               X=rbind(xx2,xx3,xx5,xx6)[,-3],R=R, 
                               w_pre=para.dd2[t-1,1:K],
                               beta0_pre=para.dd2[t-1,ind1], 
                               beta_pre=para.dd2[t-1,ind3],
                               sigma_pre=para.dd2[t-1,ind2],
                               alpha_pre=para.dd2[t-1,ind4], 
                               alpha_beta0_pre=para.dd2[t-1,ind5], 
                               alpha_sigma_pre=para.dd2[t-1,ind6],
                               mu_beta0_pre=para.dd2[t-1,ind7], 
                               sigma_beta0_pre=para.dd2[t-1,ind8], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=2, 
                               cov1=cov.d12,cov2=cov2.d12,gamma=para.dd2[t-1,(ind9-n2-n3-n5-n6)], 
                               Z=para.dd2[t-1,(ind9-n2-n3-n5-n6+1):ind9])
    
    Zm2 <- para.dd2[t,(ind9-n2-n3-n5-n6+1):ind9][(n2+n3+1):(n2+n3+n5+n6)]
    h[(n2+n3+1):(n2+n3+n5+n6),2] <- qnorm(pmin(1-0.1^10,pmax(0.1^10, ptruncnorm(rbind(DD5,DD6)[,2],a=0,b=1,
                                                                                para.dd2[t,ind1][Zm2]+
                                                                                  rbind(xx5,xx6)[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                                                                sqrt(para.dd2[t,ind2])[Zm2]))))
    
    
    
    if(t > 200){
      cov.d13 <- cov(para.dd3[3:(t-1),(K+1):(2*K)])
      cov2.d13 <- cov(para.dd3[3:(t-1),ind13])
    }
    
    h <- ifelse(is.nan(h), -3, h)
    para.dd3[t,] <- metropolis(h=h, Y=rbind(DD2,DD3,DD5,DD6)[,3], X=rbind(xx2,xx3,xx5,xx6),R=R, 
                               w_pre=para.dd3[t-1,1:K],
                               beta0_pre=para.dd3[t-1,ind1], 
                               beta_pre=para.dd3[t-1,ind13],
                               sigma_pre=para.dd3[t-1,ind12],
                               alpha_pre=para.dd3[t-1,ind14], 
                               alpha_beta0_pre=para.dd3[t-1,ind15], 
                               alpha_sigma_pre=para.dd3[t-1,ind16],
                               mu_beta0_pre=para.dd3[t-1,ind17], 
                               sigma_beta0_pre=para.dd3[t-1,ind18], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=3, 
                               cov1=cov.d13,cov2=cov2.d13,gamma=para.dd3[t-1,(ind19-n2-n3-n5-n6)], 
                               Z=para.dd3[t-1,(ind19-n2-n3-n5-n6+1):ind19])
    
    Zm3 <- para.dd3[t,(ind19-n2-n3-n5-n6+1):ind19][c(1:n2,(n2+n3+1):(n2+n3+n5))]
    h[c(1:n2,(n2+n3+1):(n2+n3+n5)),3] <- qnorm(pmin(1-0.1^10,pmax(0.1^10, ptruncnorm(rbind(DD2,DD5)[,3],a=0,b=1,
                                                                                     para.dd3[t,ind1][Zm3]+
                                                                                       rbind(xx2,xx5)%*%as.matrix(para.dd3[t,ind13]),
                                                                                     sqrt(para.dd3[t,ind12])[Zm3]))))
    if(t > 200){
      cov.d14 <- cov(para.dd4[3:(t-1),(K+1):(2*K)])
      cov2.d14 <- cov(para.dd4[3:(t-1),ind13])
    }
    
    h <- ifelse(is.nan(h), -3, h)
    para.dd4[t,] <- metropolis(h=h, Y=rbind(DD2,DD3,DD5,DD6)[,4], 
                               X=rbind(xx2,xx3,xx5,xx6),R=R, 
                               w_pre=para.dd4[t-1,1:K],
                               beta0_pre=para.dd4[t-1,ind1], 
                               beta_pre=para.dd4[t-1,ind13],
                               sigma_pre=para.dd4[t-1,ind12],
                               alpha_pre=para.dd4[t-1,ind14], 
                               alpha_beta0_pre=para.dd4[t-1,ind15], 
                               alpha_sigma_pre=para.dd4[t-1,ind16],
                               mu_beta0_pre=para.dd4[t-1,ind17], 
                               sigma_beta0_pre=para.dd4[t-1,ind18], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=4, 
                               cov1=cov.d14,cov2=cov2.d14,gamma=para.dd4[t-1,(ind19-n2-n3-n5-n6)], 
                               Z=para.dd4[t-1,(ind19-n2-n3-n5-n6+1):ind19])
    
    Zm4 <- para.dd4[t,(ind19-n2-n3-n5-n6+1):ind19][c((n2+1):(n2+n3),(n2+n3+n5+1):(n2+n3+n5+n6))]
    h[c((n2+1):(n2+n3),(n2+n3+n5+1):(n2+n3+n5+n6)),4] <- qnorm(pmin(1-0.1^10,pmax(0.1^5, ptruncnorm(rbind(DD3,DD6)[,2],a=0,b=1,
                                                                                                    para.dd4[t,ind1][Zm4]+rbind(xx3,xx6)%*%as.matrix(para.dd4[t,ind13]),
                                                                                                    sqrt(para.dd4[t,ind12])[Zm4]))))
    h <- ifelse(is.nan(h), -3, h)
    para.C[t,] <- metropolisC1(h=h, rho=para.C[t-1,7:12])
    prop1 <- para.C[t,7:12]
    
    R <- matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],
                  1,prop1[6],prop1[3],prop1[5],prop1[6],1),4,4,byrow=TRUE)
    
    
    
    h2_prop <- h2
    DD2_prop <- DD2
    
    DD2_prop[,2] <- rtruncnorm(n2,DD2[,2],sd(DD2[,2]),a=0,b=1)
    clus3 <- para.dd2[t,(ind9-n2-n3-n5-n6+1):ind9][1:n2]
    clusy <- para.y1[t,(indy91-n2-n5+1):indy91][1:n2]
    h2_prop[,2] <-qnorm(pmin(1-0.1^10,pmax(0.1^10,
                                           ptruncnorm(DD2_prop[,2],a=0,b=1,
                                                      para.dd2[t,ind1][clus3]+
                                                        xx2[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                                      sqrt(para.dd2[t,ind2])[clus3]))))
    
    
    h2_prop <- ifelse(is.nan(h2_prop), -3, h2_prop)
    
    
    rat1 <- dmnorm(h2_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2_prop[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus3]+
                                    xx2[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus3])))+
      log(pmax(0.1^100,dnorm(y2[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx2,DD2_prop,aa2)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2[c,2],DD2_prop[c,2], 
                                              sd(DD2[,2]),a=0,b=1)))
    
    rat2 <- dmnorm(h2,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus3]+
                                    xx2[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus3])))+
      log(pmax(0.1^100,dnorm(y2[,"y"],para.y2[t,indy1][clusy]+
                               cbind(xx2,DD2,aa2)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy])))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2_prop[c,2],DD2[c,2], 
                                              sd(DD2[,2]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc2 <- rep(0,n2)
    
    for(ii in 1:n2){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD2_prop[ii,2] <- DD2[ii,2]
        h2_prop[ii,2] <- h2[ii,2]
      }else{
        DD2[ii,2] <- DD2_prop[ii,2]
        h2[ii,2] <- h2_prop[ii,2]
        acc2[ii] <- 1
      }
    }
    
    DD2_prop[,4] <- rtruncnorm(n2,DD2[,4],sd(DD2[,4]),a=0,b=1)
    clus03 <- para.dd4[t,(ind19-n2-n3-n5-n6+1):ind19][1:n2]
    clusy <- para.y1[t,(indy91-n2-n5+1):indy91][1:n2]
    h2_prop[,4] <-qnorm(pmin(1-0.1^10,pmax(0.1^10,
                                           ptruncnorm(DD2_prop[,4],a=0,b=1,
                                                      para.dd4[t,ind1][clus03]+
                                                        xx2%*%as.matrix(para.dd4[t,ind13]),
                                                      sqrt(para.dd4[t,ind12])[clus03]))))
    
    h2_prop <- ifelse(is.nan(h2_prop), -3, h2_prop)
    
    
    rat1 <- dmnorm(h2_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2_prop[,4],a=0,b=1,
                                  para.dd4[t,ind1][clus03]+
                                    xx2%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus03])))+
      log(pmax(0.1^100,dnorm(y2[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx2,DD2_prop,aa2)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2[c,4],DD2_prop[c,4], 
                                              sd(DD2[,4]),a=0,b=1)))
    
    rat2 <- dmnorm(h2,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2[,4],a=0,b=1,
                                  para.dd4[t,ind1][clus03]+
                                    xx2%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus03])))+
      log(pmax(0.1^100,dnorm(y2[,"y"],para.y2[t,indy1][clusy]+
                               cbind(xx2,DD2,aa2)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy])))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2_prop[c,4],DD2[c,4], 
                                              sd(DD2[,4]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc2 <- rep(0,n2)
    
    for(ii in 1:n2){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD2_prop[ii,4] <- DD2[ii,4]
        h2_prop[ii,4] <- h2[ii,4]
      }else{
        DD2[ii,4] <- DD2_prop[ii,4]
        h2[ii,4] <- h2_prop[ii,4]
        acc2[ii] <- 1
      }
    }
    
    h3_prop <- h3
    DD3_prop <- DD3
    
    DD3_prop[,2] <- rtruncnorm(n3,DD3[,2],sd(DD3[,2]),a=0,b=1)
    clus4 <- para.dd2[t,(ind9-n2-n3-n5-n6+1):ind9][(n2+1):(n2+n3)]
    clusy2 <- para.y2[t,(indy92-n3-n6+1):indy92][1:n3]
    h3_prop[,2] <-qnorm(pmin(1-0.1^10,pmax(0.1^10,
                                           ptruncnorm(DD3_prop[,2],a=0,b=1,
                                                      para.dd2[t,ind1][clus4]+
                                                        xx3[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                                      sqrt(para.dd2[t,ind2])[clus4]))))
    
    
    #h3_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #pnorm(DD3_prop[,3],para.dd3[t,ind1][clus5],
    #sqrt(para.dd3[t,ind2])[clus5]))))
    
    h3_prop <- ifelse(is.nan(h3_prop), -3, h3_prop)
    
    rat1 <- dmnorm(h3_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3_prop[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus4]+
                                    xx3[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus4])))+
      log(pmax(0.1^100,dnorm(y3[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx3,DD3_prop,aa3)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n3, function(c) log(dtruncnorm(DD3[c,2],DD3_prop[c,2], 
                                              sd(DD3[,2]),a=0,b=1)))
    
    rat2 <- dmnorm(h3,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus4]+
                                    xx3[,-3]%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus4])))+
      log(pmax(0.1^100,dnorm(y3[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx3,DD3,aa3)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n3, function(c) log(dtruncnorm(DD3_prop[c,2],DD3[c,2], 
                                              sd(DD3[,2]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc3 <- rep(0,n3)
    
    for(ii in 1:n3){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD3_prop[ii,2] <- DD3[ii,2]
        h3_prop[ii,2] <- h3[ii,2]
      }else{
        DD3[ii,2] <- DD3_prop[ii,2]
        h3[ii,2] <- h3_prop[ii,2]
        acc3[ii] <- 1
      }
    }
    
    DD3_prop[,3] <- rtruncnorm(n3,DD3[,3],sd(DD3[,3]),a=0,b=1)
    clus5 <- para.dd3[t,(ind19-n2-n3-n5-n6+1):ind19][(n2+1):(n2+n3)]
    h3_prop[,3] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD3_prop[,3],a=0,b=1,
                                                     para.dd3[t,ind1][clus5]+
                                                       xx3%*%as.matrix(para.dd3[t,ind13]),
                                                     sqrt(para.dd3[t,ind12])[clus5]))))
    
    h3_prop <- ifelse(is.nan(h3_prop), -3, h3_prop)
    
    rat1 <- dmnorm(h3_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3_prop[,3],a=0,b=1,
                                  para.dd3[t,ind1][clus5]+
                                    xx3%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus5])))+
      log(pmax(0.1^100,dnorm(y3[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx3,DD3_prop,aa3)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n3, function(c) log(dtruncnorm(DD3[c,3],DD3_prop[c,3], 
                                              sd(DD3[,3]),a=0,b=1)))
    
    rat2 <- dmnorm(h3,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3[,3],a=0,b=1,
                                  para.dd3[t,ind1][clus5]+
                                    xx3%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus5])))+
      log(pmax(0.1^100,dnorm(y3[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx3,DD3,aa3)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n3, function(c) log(dtruncnorm(DD3_prop[c,3],DD3[c,3], 
                                              sd(DD3[,3]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc3 <- rep(0,n3)
    
    for(ii in 1:n3){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD3_prop[ii,3] <- DD3[ii,3]
        h3_prop[ii,3] <- h3[ii,3]
      }else{
        DD3[ii,3] <- DD3_prop[ii,3]
        h3[ii,3] <- h3_prop[ii,3]
        acc3[ii] <- 1
      }
    }
    
    
    h5_prop <- h5
    DD5_prop <- DD5
    
    DD5_prop[,1] <- rtruncnorm(n5,DD5[,1],sd(DD5[,1]),a=0,b=1)
    clus8 <- para.dd1[t,(ind9-n2-n3-n5-n6+1):ind9][(n2+n3+1):(n2+n3+n5)]
    clusy <- para.y1[t,(indy91-n2-n5+1):indy91][(n2+1):(n2+n5)]
    h5_prop[,1] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD5_prop[,1],a=0,b=1,para.dd1[t,ind1][clus8]+
                                                       xx5[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                                     sqrt(para.dd1[t,ind2])[clus8]))))
    
    
    h5_prop <- ifelse(is.nan(h5_prop), -3, h5_prop)
    
    
    rat1 <- dmnorm(h5_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5_prop[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus8]+
                                    xx5[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus8])))+
      log(pmax(0.1^100,dnorm(y5[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx5,DD5_prop,aa5)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5[c,1],DD5_prop[c,1], 
                                              sd(DD5[,1]),a=0,b=1)))
    
    rat2 <- dmnorm(h5,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus8]+
                                    xx5[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus8])))+
      log(pmax(0.1^100,dnorm(y5[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx5,DD5,aa5)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5_prop[c,1],DD5[c,1], 
                                              sd(DD5[,1]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc5 <- rep(0,n5)
    
    for(ii in 1:n5){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD5_prop[ii,1] <- DD5[ii,1]
        h5_prop[ii,1] <- h5[ii,1]
      }else{
        DD5[ii,1] <- DD5_prop[ii,1]
        h5[ii,1] <- h5_prop[ii,1]
        acc5[ii] <- 1
      }
    }
    
    DD5_prop[,4] <- rtruncnorm(n5,DD5[,4],sd(DD5[,4]),a=0,b=1)
    clus08 <- para.dd4[t,(ind19-n2-n3-n5-n6+1):ind19][(n2+n3+1):(n2+n3+n5)]
    clusy <- para.y1[t,(indy91-n2-n5+1):indy91][(n2+1):(n2+n5)]
    h5_prop[,4] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
                                          ptruncnorm(DD5_prop[,4],para.dd4[t,ind1][clus08]+
                                                       xx5%*%as.matrix(para.dd4[t,ind13]),
                                                     sqrt(para.dd4[t,ind12])[clus08]))))
    
    h5_prop <- ifelse(is.nan(h5_prop), -3, h5_prop)
    
    
    rat1 <- dmnorm(h5_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5_prop[,4],a=0,b=1,
                                  para.dd4[t,ind1][clus08]+
                                    xx5%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus08])))+
      log(pmax(0.1^100,dnorm(y5[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx5,DD5_prop,aa5)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5[c,4],DD5_prop[c,4], 
                                              sd(DD5[,4]),a=0,b=1)))
    
    rat2 <- dmnorm(h5,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5[,4],a=0,b=1,
                                  para.dd4[t,ind1][clus08]+
                                    xx5%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus08])))+
      log(pmax(0.1^100,dnorm(y5[,"y"],para.y1[t,indy1][clusy]+
                               cbind(xx5,DD5,aa5)%*%as.matrix(para.y1[t,indy3]),
                             sqrt(para.y1[t,indy2])[clusy])))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5_prop[c,4],DD5[c,4], 
                                              sd(DD5[,4]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc5 <- rep(0,n5)
    
    for(ii in 1:n5){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD5_prop[ii,4] <- DD5[ii,4]
        h5_prop[ii,4] <- h5[ii,4]
      }else{
        DD5[ii,4] <- DD5_prop[ii,4]
        h5[ii,4] <- h5_prop[ii,4]
        acc5[ii] <- 1
      }
    }
    
    h6_prop <- h6
    DD6_prop <- DD6
    
    DD6_prop[,1] <- rtruncnorm(n6,DD6[,1],sd(DD6[,1]),a=0,b=1)
    clus9 <- para.dd1[t,(ind9-n2-n3-n5-n6+1):ind9][(n2+n3+n5+1):(n2+n3+n5+n6)]
    clusy2 <- para.y2[t,(indy92-n3-n6+1):indy92][(n3+1):(n3+n6)]
    h6_prop[,1] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD6_prop[,1],a=0,b=1,
                                                     para.dd1[t,ind1][clus9]+
                                                       xx6[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                                     sqrt(para.dd1[t,ind2])[clus9]))))
    
    
    h6_prop <- ifelse(is.nan(h6_prop), -3, h6_prop)
    
    rat1 <- dmnorm(h6_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6_prop[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus9]+
                                    xx6[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus9])))+
      log(pmax(0.1^100,dnorm(y6[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx6,DD6_prop,aa6)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n6, function(c) log(dtruncnorm(DD6[c,1],DD6_prop[c,1], 
                                              sd(DD6[,1]),a=0,b=1)))
    
    rat2 <- dmnorm(h6,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6[,1],a=0,b=1,para.dd1[t,ind1][clus9]+
                                    xx6[,-3]%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus9])))+
      log(pmax(0.1^100,dnorm(y6[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx6,DD6,aa6)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n6, function(c) log(dtruncnorm(DD6_prop[c,1],DD6[c,1], 
                                              sd(DD6[,1]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc6 <- rep(0,n6)
    
    for(ii in 1:n6){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD6_prop[ii,1] <- DD6[ii,1]
        h6_prop[ii,1] <- h6[ii,1]
      }else{
        DD6[ii,1] <- DD6_prop[ii,1]
        h6[ii,1] <- h6_prop[ii,1]
        acc6[ii] <- 1
      }
    }
    
    h6_prop <- h6
    DD6_prop <- DD6
    
    DD6_prop[,3] <- rtruncnorm(n6,DD6[,3],sd(DD6[,3]),a=0,b=1)
    clusy2 <- para.y2[t,(indy92-n3-n6+1):indy92][(n3+1):(n3+n6)]
    
    clus09 <- para.dd3[t,(ind19-n2-n3-n5-n6+1):ind19][(n2+n3+n5+1):(n2+n3+n5+n6)]
    h6_prop[,3] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD6_prop[,3],a=0,b=1,
                                                     para.dd3[t,ind1][clus09]+
                                                       xx6%*%as.matrix(para.dd3[t,ind13]),
                                                     sqrt(para.dd3[t,ind12])[clus09]))))
    
    h6_prop <- ifelse(is.nan(h6_prop), -3, h6_prop)
    
    rat1 <- dmnorm(h6_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6_prop[,3],a=0,b=1,
                                  para.dd3[t,ind1][clus09]+
                                    xx6%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus09])))+
      log(pmax(0.1^100,dnorm(y6[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx6,DD6_prop,aa6)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n6, function(c) log(dtruncnorm(DD6[c,3],DD6_prop[c,3], 
                                              sd(DD6[,3]),a=0,b=1)))
    
    rat2 <- dmnorm(h6,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6[,3],a=0,b=1,para.dd3[t,ind1][clus09]+
                                    xx6%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus09])))+
      log(pmax(0.1^100,dnorm(y6[,"y"],para.y2[t,indy1][clusy2]+
                               cbind(xx6,DD6,aa6)%*%as.matrix(para.y2[t,indy3]),
                             sqrt(para.y2[t,indy2])[clusy2])))+
      sapply(1:n6, function(c) log(dtruncnorm(DD6_prop[c,3],DD6[c,3], 
                                              sd(DD6[,3]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc6 <- rep(0,n6)
    
    for(ii in 1:n6){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD6_prop[ii,3] <- DD6[ii,3]
        h6_prop[ii,3] <- h6[ii,3]
      }else{
        DD6[ii,3] <- DD6_prop[ii,3]
        h6[ii,3] <- h6_prop[ii,3]
        acc6[ii] <- 1
      }
    }
    h <- rbind(h2,h3,h5,h6)
    D <- rbind(DD2,DD3,DD5,DD6)
    ME[[t]] <- list(D=D,h=h,R=R)
    yy_2 <- yy_3 <- yy_5 <- yy_6 <- c()
    for(l in 1:(n2+n3)){
      yy_2[l] <- sum(para.y1[t,1:K1]*(para.y1[t,indy1]+sum(c(x[n1+l,],D[l,],1)* para.y1[t,indy3])))
      
      
      yy_3[l]<- sum(para.y2[t,1:K1]*(para.y2[t,indy1]+sum(c(x[n1+l,],D[l,],1)* para.y2[t,indy3])))
    }
    for(l in 1:(n5+n6)){
      yy_5[l] <- sum(para.y1[t,1:K1]*(para.y1[t,indy1]+sum(c(x[n1+n2+n3+n4+l,],D[n2+n3+l,],-1)* para.y1[t,indy3])))
      
      
      yy_6[l]<- sum(para.y2[t,1:K1]*(para.y2[t,indy1]+sum(c(x[n1+n2+n3+n4+l,],D[n2+n3+l,],-1)* para.y2[t,indy3])))
    }
    
    yyi1 <- ifelse(yy_2>yy_3,1,0)
    yyi2 <- ifelse(yy_5>yy_6,1,0)
    x23 <- rbind(xx2,xx3)
    psu1 <- psu2 <- c()
    for(l in 1:length(yyi1)){
      if(yyi1[l]==1){
        psu1[l]<- sum(para.y1[t,1:K]*(para.y1[t,indy1]+sum(c(x[n1+l,],D[l,],1)* para.y1[t,indy3])))
      }
      else{
        psu1[l]<- sum(para.y2[t,1:K]*(para.y2[t,indy1]+sum(c(x[n1+l,],D[l,],1)* para.y2[t,indy3])))
        
      }
    }
    psu[(n1+1):(n1+n2+n3)] <- psu1
    x56 <- rbind(xx5,xx6)
    for(l in 1:length(yyi2)){
      if(yyi2[l]==1){
        psu2[l]<- sum(para.y1[t,1:K]*(para.y1[t,indy1]+sum(c(x[n1+n2+n3+n4+l,],D[n2+n3+l,],-1)* para.y1[t,indy3])))
      }
      else{
        psu2[l]<- sum(para.y2[t,1:K]*(para.y2[t,indy1]+sum(c(x[n1+n2+n3+n4+l,],D[n2+n3+l,],-1)* para.y2[t,indy3])))
        
      }
    }
    psu[(n1+n2+n3+n4+1):n] <- psu2
    DD <- rbind(DD1,D[1:(n2+n3),],DD4,D[(n2+n3+1):(n2+n3+n5+n6),])
    if(t==3){
      R1 <- diag(2)
      para.d1 <- matrix(nrow=mcmc,ncol=K+K+2+K+5+1+n)
      para.d1[1,] <-para.d1[2,] <- c(rep(1/K,K),rep(coefficients(lm(DD[1:(n1+n2+n3),1]~x[1:(n1+n2+n3),1]+
                                                                      x[1:(n1+n2+n3),2]))[1],K),
                                     coefficients(lm(DD[1:(n1+n2+n3),1]~x[1:(n1+n2+n3),1]+
                                                       x[1:(n1+n2+n3),2]))[-1],
                                     rep(1,K),1,1,1,mean(DD[1:(n1+n2+n3),1]),1,1, 
                                     sample(1:K,n, replace=TRUE))
      
      para.d2 <- matrix(nrow=mcmc,ncol=K+K+2+K+5+1+n)
      para.d2[1,] <-para.d2[2,] <- c(rep(1/K,K),rep(coefficients(lm(DD[(n1+n2+n3+1):n,1]~x[(n1+n2+n3+1):n,1]+
                                                                      x[(n1+n2+n3+1):n,2]))[1],K),
                                     coefficients(lm(DD[(n1+n2+n3+1):n,1]~x[(n1+n2+n3+1):n,1]+
                                                       x[(n1+n2+n3+1):n,2]))[-1],
                                     rep(1,K),1,1,1,mean(DD[(n1+n2+n3+1):n,2]),1,1, 
                                     sample(1:K,n, replace=TRUE))
      
      
      
      para.C1<-matrix(nrow = mcmc, ncol = 2)
      para.C1[1,] <- para.C1[2,] <- c(NA,0)
      
      indd1 <- (K+1):(2*K) # intercept
      indd3<- (2*K+1):(2*K+2)
      indd2 <- (2*K+3):(3*K+2) # variance
      indd4 <- 3*K+3 # alpha
      indd5 <- indd4+1 # alpha_beta0
      indd6 <- indd5+1 # alpha_sigma
      indd7 <- indd6+1 # mu_beta0
      indd8 <- indd7+1 # sigma_beta0
      
      indd9 <- length(para.d1[1,])
      
      
      h11 <- DD[1:(n1+n2+n3),]; h12 <- DD[(n1+n2+n3+1):n,]
      
      ME1 <- list()
      par.y1 <- matrix(nrow=mcmc,ncol=K1+K1+4+K1+5+1+n1+n2+n3)
      par.y1[1,] <-par.y1[2,] <- c(rep(1/K1,K1),rep(coefficients(lm(psu[1:(n1+n2+n3)]~x[1:(n1+n2+n3),1]+
                                                                      x[1:(n1+n2+n3),2]+DD[1:(n1+n2+n3),1]+DD[1:(n1+n2+n3),2]))[1],K1),
                                   coefficients(lm(psu[1:(n1+n2+n3)]~x[1:(n1+n2+n3),1]+
                                                     x[1:(n1+n2+n3),2]+DD[1:(n1+n2+n3),1]+DD[1:(n1+n2+n3),2]))[-1],
                                   rep(1,K1),1,1,1,mean(psu[1:(n1+n2+n3)]),1,1, 
                                   sample(1:K1,n1+n2+n3, replace=TRUE))
      
      par.y2 <- matrix(nrow=mcmc,ncol=K1+K1+4+K1+5+1+n4+n5+n6)
      par.y2[1,] <-par.y2[2,] <- c(rep(1/K1,K1),rep(coefficients(lm(psu[(n1+n2+n3+1):n]~x[(n1+n2+n3+1):n,1]+
                                                                      x[(n1+n2+n3+1):n,2]+DD[(n1+n2+n3+1):n,1]+DD[(n1+n2+n3+1):n,2]))[1],K1),
                                   coefficients(lm(psu[(n1+n2+n3+1):n]~x[(n1+n2+n3+1):n,1]+
                                                     x[(n1+n2+n3+1):n,2]+DD[(n1+n2+n3+1):n,1]+DD[(n1+n2+n3+1):n,2]))[-1],
                                   rep(1,K1),1,1,1,mean(psu[(n1+n2+n3+1):n]),1,1, 
                                   sample(1:K1,n4+n5+n6, replace=TRUE))
      
      
      indyy1 <- (K1+1):(2*K1) # intercept
      indyy3<- (2*K1+1):(2*K1+4)
      indyy2 <- (2*K1+5):(3*K1+4) # variance
      indyy4 <- 3*K1+5 # alpha
      indyy5 <- indyy4+1 # alpha_beta0
      indyy6 <- indyy5+1 # alpha_sigma
      indyy7 <- indyy6+1 # mu_beta0
      indyy8 <- indyy7+1 # sigma_beta0
      
      indyy91 <- length(par.y1[1,]); indyy92 <- length(par.y2[1,])
    }
    
    cov.d1 <- diag(vcov(lm(DD[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                             rbind(xx1,xx2,xx3)[,2]))[1,1],K)
    cov.d2 <- diag(vcov(lm(DD[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                             rbind(xx4,xx5,xx6)[,2]))[1,1],K)
    
    cov2.d1 <- diag(diag(vcov(lm(DD[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                                   rbind(xx1,xx2,xx3)[,2]))[-1,-1]),2)
    cov2.d2 <- diag(diag(vcov(lm(DD[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                                   rbind(xx4,xx5,xx6)[,2]))[-1,-1]),2)
    
    co.y1 <- diag(vcov(lm(psu[1:(n1+n2+n3)]~x[1:(n1+n2+n3),1]+
                            x[1:(n1+n2+n3),2]+DD[1:(n1+n2+n3),1]+DD[1:(n1+n2+n3),2]))[1,1],K1)
    
    co2.y1 <- diag(diag(vcov(lm(psu[1:(n1+n2+n3)]~x[1:(n1+n2+n3),1]+
                                  x[1:(n1+n2+n3),2]+DD[1:(n1+n2+n3),1]+DD[1:(n1+n2+n3),2]))[-1,-1]),4)
    
    co.y2 <- diag(vcov(lm(psu[(n1+n2+n3+1):n]~x[(n1+n2+n3+1):n,1]+
                            x[(n1+n2+n3+1):n,2]+DD[(n1+n2+n3+1):n,1]+DD[(n1+n2+n3+1):n,2]))[1,1],K1)
    
    co2.y2 <- diag(diag(vcov(lm(psu[(n1+n2+n3+1):n]~x[(n1+n2+n3+1):n,1]+
                                  x[(n1+n2+n3+1):n,2]+DD[(n1+n2+n3+1):n,1]+DD[(n1+n2+n3+1):n,2]))[-1,-1]),4)
    acc2 <-acc3 <-acc5<-acc6<-0
    
    if(t > 200){
      co.y1 <- cov(par.y1[3:(t-1),(K1+1):(2*K1)])
      co2.y1 <- cov(par.y1[3:(t-1),indyy3])
      co.y2 <- cov(par.y2[3:(t-1),(K1+1):(2*K1)])
      co2.y2 <- cov(par.y2[3:(t-1),indyy3])
    }
    par.y1[t,] <- metropolisy(Y=psu[1:(n1+n2+n3)], X=cbind(x[1:(n1+n2+n3),1:2],DD[1:(n1+n2+n3),1:2]),
                              w_pre=par.y1[t-1,1:K1],
                              beta0_pre=par.y1[t-1,indyy1], 
                              beta_pre=par.y1[t-1,indyy3],
                              sigma_pre=par.y1[t-1,indyy2],
                              alpha_pre=par.y1[t-1,indyy4], 
                              alpha_beta0_pre=par.y1[t-1,indyy5], 
                              alpha_sigma_pre=par.y1[t-1,indyy6],
                              mu_beta0_pre=par.y1[t-1,indyy7], 
                              sigma_beta0_pre=par.y1[t-1,indyy8], K=K1, 
                              eps=0.1,del1=15, del2=10, del3=20, 
                              cov1=co.y1,cov2=co2.y1,gamma=par.y1[t-1,(indyy91-n1-n2-n3)], 
                              Z=par.y1[t-1,(indyy91-n1-n2-n3+1):indyy91])
    
    par.y2[t,] <- metropolisy(Y=psu[(n1+n2+n3+1):n], X=cbind(x[(n1+n2+n3+1):n,1:2],DD[(n1+n2+n3+1):n,1:2]),
                              w_pre=par.y2[t-1,1:K1],
                              beta0_pre=par.y2[t-1,indyy1], 
                              beta_pre=par.y2[t-1,indyy3],
                              sigma_pre=par.y2[t-1,indyy2],
                              alpha_pre=par.y2[t-1,indyy4], 
                              alpha_beta0_pre=par.y2[t-1,indyy5], 
                              alpha_sigma_pre=par.y2[t-1,indyy6],
                              mu_beta0_pre=par.y2[t-1,indyy7], 
                              sigma_beta0_pre=par.y2[t-1,indyy8], K=K1, 
                              eps=0.1,del1=15, del2=10, del3=20, 
                              cov1=co.y2,cov2=co2.y2,gamma=par.y2[t-1,(indyy92-n4-n5-n6)], 
                              Z=par.y2[t-1,(indyy92-n4-n5-n6+1):indyy92])
    
    if(t > 200){
      cov.d1 <- cov(para.d1[3:(t-1),(K+1):(2*K)])
      cov2.d1 <- cov(para.d1[3:(t-1),ind3])
    }
    
    hh <- rbind(h11,h12)
    hh <- ifelse(is.nan(hh), -3, hh)
    para.d1[t,] <- metropolis(h=hh[,1:2], Y=DD[,1], X=x[,1:2],
                              R=R1, w_pre=para.d1[t-1,1:K],
                              beta0_pre=para.d1[t-1,indd1], 
                              beta_pre=para.d1[t-1,indd3],
                              sigma_pre=para.d1[t-1,indd2],
                              alpha_pre=para.d1[t-1,indd4], 
                              alpha_beta0_pre=para.d1[t-1,indd5], 
                              alpha_sigma_pre=para.d1[t-1,indd6],
                              mu_beta0_pre=para.d1[t-1,indd7], 
                              sigma_beta0_pre=para.d1[t-1,indd8], K=K, 
                              eps=0.1,del1=15, del2=10, del3=20, index=1, 
                              cov1=cov.d1,cov2=cov2.d1,gamma=para.d1[t-1,(indd9-n)], 
                              Z=para.d1[t-1,(indd9-n+1):indd9])
    
    Zm1 <- para.d1[t,(indd9-n+1):indd9][1:(n1+n2+n3)]
    hh[1:(n1+n2+n3),1] <- qnorm(pmin(1-0.1^10,pmax(0.1^10, ptruncnorm(DD[1:(n1+n2+n3),1],a=0,b=1,
                                                                      para.d1[t,indd1][Zm1]+
                                                                        x[1:(n1+n2+n3),-3]%*%as.matrix(para.d1[t,indd3]),
                                                                      sqrt(para.d1[t,indd2])[Zm1]))))
    
    
    
    if(t > 200){
      cov.d12 <- cov(para.d2[3:(t-1),(K+1):(2*K)])
      cov2.d12 <- cov(para.d2[3:(t-1),ind3])
    }
    
    hh <- ifelse(is.nan(hh), -3, hh)
    para.d2[t,] <- metropolis(h=hh[,1:2], Y=DD[,2], X=x[,-3],R=R1, 
                              w_pre=para.d2[t-1,1:K],
                              beta0_pre=para.d2[t-1,indd1], 
                              beta_pre=para.d2[t-1,indd3],
                              sigma_pre=para.d2[t-1,indd2],
                              alpha_pre=para.d2[t-1,indd4], 
                              alpha_beta0_pre=para.d2[t-1,indd5], 
                              alpha_sigma_pre=para.d2[t-1,indd6],
                              mu_beta0_pre=para.d2[t-1,indd7], 
                              sigma_beta0_pre=para.d2[t-1,indd8], K=K, 
                              eps=0.1,del1=15, del2=10, del3=20, index=2, 
                              cov1=cov.d12,cov2=cov2.d12,gamma=para.d2[t-1,(indd9-n)], 
                              Z=para.d2[t-1,(indd9-n+1):indd9])
    
    Zm2 <- para.d2[t,(indd9-n+1):indd9][(n1+n2+n3+1):n]
    hh[(n1+n2+n3+1):n,2] <- qnorm(pmin(1-0.1^10,pmax(0.1^10, ptruncnorm(DD[(n1+n2+n3+1):n,2],a=0,b=1,
                                                                        para.d2[t,indd1][Zm2]+
                                                                          x[(n1+n2+n3+1):n,-3]%*%as.matrix(para.d2[t,indd3]),
                                                                        sqrt(para.d2[t,indd2])[Zm2]))))
    
    
    
    hh <- ifelse(is.nan(hh), -3, hh)
    para.C1[t,] <- metropolisC2(h=hh[,1:2], rho=para.C1[t-1,2])
    prop1 <- para.C1[t,2]
    
    R1 <- matrix(c(1,prop1,prop1,1),2,2,byrow=TRUE)
    
    
    
    
    h11_prop <- h11
    DD1_prop <- DD1
    
    DD1_prop[,2] <- rtruncnorm(n1,DD1[,2],sd(DD1[,2]),a=0,b=1)
    clus1 <- para.d2[t,(indd9-n+1):indd9][1:n1]
    clusy1 <- par.y1[t,(indyy91-n1-n2-n3+1):indyy91][1:n1]
    h11_prop[1:n1,2] <-qnorm(pmin(1-0.1^10,pmax(0.1^10,
                                                ptruncnorm(DD1_prop[,2],a=0,b=1,
                                                           para.d2[t,indd1][clus1]+
                                                             xx1[,-3]%*%as.matrix(para.d2[t,indd3]),
                                                           sqrt(para.d2[t,indd2])[clus1]))))
    
    h11_prop <- ifelse(is.nan(h11_prop), -3, h11_prop)
    #clus2 <- para.dd3[t,(ind9-n+1):ind9][1:n1]
    #h1_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #                                     pnorm(DD1_prop[,3],para.dd3[t,ind1][clus2],
    #                                           sqrt(para.dd3[t,ind2])[clus2]))))
    
    
    rat1 <- dmnorm(h11_prop[1:n1,1:2],rep(0,2), R1, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD1_prop[,2],a=0,b=1,para.d2[t,indd1][clus1]+
                                    xx1[,-3]%*%as.matrix(para.d2[t,indd3]),
                                  sqrt(para.d2[t,indd2])[clus1])))+
      #log(pmax(0.1^100,dnorm(DD1_prop[,3])))+
      log(pmax(0.1^100,dnorm(psu[1:n1],par.y1[t,indyy1][clusy1]+
                               cbind(xx1[,1:2],DD1_prop[,1:2])%*%as.matrix(par.y1[t,indyy3]),
                             sqrt(par.y1[t,indyy2])[clusy1])))+
      sapply(1:n1, function(c) log(dtruncnorm(DD1[c,2],DD1_prop[c,2], 
                                              sd(DD1[,2]),a=0,b=1)))
    
    rat2 <- dmnorm(h11[1:n1,1:2],rep(0,2), R1, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD1[,2],a=0,b=1,para.d2[t,indd1][clus1]+
                                    xx1[,-3]%*%as.matrix(para.d2[t,indd3]),
                                  sqrt(para.d2[t,indd2])[clus1])))+
      #log(pmax(0.1^100,dnorm(DD1[,3])))+
      log(pmax(0.1^100,dnorm(psu[1:n1],par.y1[t,indyy1][clusy1]+
                               cbind(xx1[,1:2],DD1[,1:2])%*%as.matrix(par.y1[t,indyy3]),
                             sqrt(par.y1[t,indyy2])[clusy1])))+
      sapply(1:n1, function(c) log(dtruncnorm(DD1_prop[c,2],DD1[c,2], 
                                              sd(DD1[,2]), a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc1 <- rep(0,n1)
    
    for(ii in 1:n1){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD1_prop[ii,2] <- DD1[ii,2]
        h11_prop[ii,2] <- h11[ii,2]
      }else{
        DD1[ii,2] <- DD1_prop[ii,2]
        h11[ii,2] <- h11_prop[ii,2]
        acc1[ii] <- 1
      }
    }
    
    h12_prop<-h12; DD4_prop<-DD4
    
    DD4_prop[,1] <- rtruncnorm(n4,DD4[,1],sd(DD4[,1]),a=0,b=1)
    clus6 <- para.d1[t,(indd9-n+1):indd9][(1+n1+n2+n3):(n1+n2+n3+n4)]
    clusy2 <- par.y2[t,(indyy92-n4-n5-n6+1):indyy92][1:n4]
    h12_prop[1:n4,1] <-qnorm(pmin(1-0.1^10,pmax(0.1^10,ptruncnorm(DD4_prop[,1],a=0,b=1,
                                                                  para.d1[t,indd1][clus6]+
                                                                    xx4[,-3]%*%as.matrix(para.d1[t,indd3]),
                                                                  sqrt(para.d1[t,indd2])[clus6]))))
    
    #clus7 <- para.dd3[t,(ind9-n+1):ind9][1:n1][(1+n1+n2+n3):(n1+n2+n3+n4)]
    #h4_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #                                   pnorm(DD4_prop[,3],para.dd3[t,ind1][clus7],
    #                                        sqrt(para.dd3[t,ind2])[clus7]))))
    h12_prop <- ifelse(is.nan(h12_prop), -3, h12_prop)
    
    rat1 <- dmnorm(h12_prop[1:n4,1:2],c(0,0), R1[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD4_prop[,1],a=0,b=1,
                                  para.d1[t,indd1][clus6]+
                                    xx4[,-3]%*%as.matrix(para.d1[t,indd3]),
                                  sqrt(para.d1[t,indd2])[clus6])))+
      #log(pmax(0.1^100,dnorm(DD4_prop[,3])))+
      log(pmax(0.1^100,dnorm(psu[(n1+n2+n3+1):(n1+n2+n3+n4)],par.y2[t,indyy1][clusy2]+
                               cbind(xx4[,1:2],DD4_prop[,1:2])%*%as.matrix(par.y2[t,indyy3]),
                             sqrt(par.y2[t,indyy2])[clusy2])))+
      sapply(1:n4, function(c) log(dtruncnorm(DD4[c,1],DD4_prop[c,1], 
                                              sd(DD4[,1]),a=0,b=1)))
    
    rat2 <- dmnorm(h12[1:n4,1:2],c(0,0), R1[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD4[,1],a=0,b=1,para.d1[t,ind1][clus6]+
                                    xx4[,-3]%*%as.matrix(para.d1[t,indd3]),
                                  sqrt(para.d1[t,indd2])[clus6])))+
      #log(pmax(0.1^100,dnorm(DD4[,3])))+
      log(pmax(0.1^100,dnorm(psu[(n1+n2+n3+1):(n1+n2+n3+n4)],par.y2[t,indyy1][clusy2]+
                               cbind(xx4[,1:2],DD4[,1:2])%*%as.matrix(par.y2[t,indyy3]),
                             sqrt(par.y2[t,indyy2])[clusy2])))+
      sapply(1:n4, function(c) log(dtruncnorm(DD4_prop[c,1],DD4[c,1], 
                                              sd(DD4[,1]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc2 <- rep(0,n4)
    
    for(ii in 1:n4){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD4_prop[ii,1] <- DD4[ii,1]
        h12_prop[ii,1] <- h12[ii,1]
      }else{
        DD4[ii,1] <- DD4_prop[ii,1]
        h12[ii,1] <- h12_prop[ii,1]
        acc2[ii] <- 1
      }
    }
    DD <- rbind(DD1,DD2,DD3,DD4,DD5,DD6)
    hh <- rbind(h11,h12)
    ME1[[t]] <- list(R1=R1,DD=DD,hh=hh)
  }
  
  #Fixing the covariates at their median
  mdm <- apply(x,2,median)
  
  yy21 <- yy22 <- matrix(0,nrow(Ds),mcmc-2)
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yy21[l,j] <- sum(para.y1[j+2,1:K1]*(para.y1[j+2,indy1]+sum(c(mdm,0,0,unlist(Ds[l,]),1)* para.y1[j+2,indy3])))
      
    }
  }
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yy22[l,j] <- sum(para.y2[j+2,1:K1]*(para.y2[j+2,indy1]+sum(c(mdm,0,0,unlist(Ds[l,]),1)* para.y2[j+2,indy3])))
      
      
    }
  }
  
  
  dg21 <- matrix(0,nrow(Ds),(mcmc-2))
  for(i in 1:nrow(Ds)){
    for(j in 1:(mcmc-2)){
      dg21[i,j] <- ifelse(yy21[i,j]>yy22[i,j],1,0)
    }
  }
  
  yy23 <- yy24 <- matrix(0,nrow(Ds),(mcmc-2))
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yy23[l,j] <- sum(para.y1[j+2,1:K1]*(para.y1[j+2,indy1]+sum(c(mdm,0,0,unlist(Ds[l,]),-1)* para.y1[j+2,indy3])))
      
    }
  }
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yy24[l,j] <- sum(para.y2[j+2,1:K1]*(para.y2[j+2,indy1]+sum(c(mdm,0,0,unlist(Ds[l,]),-1)* para.y2[j+2,indy3])))
      
      
    }
  }
  
  
  dg22 <- matrix(0,nrow(Ds),(mcmc-2))
  for(i in 1:nrow(Ds)){
    for(j in 1:(mcmc-2)){
      dg22[i,j] <- ifelse(yy23[i,j]>yy24[i,j],1,0)
    }
  }
  
  yyj1 <- yyj2 <- matrix(0,nrow(Ds),mcmc-2)
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yyj1[l,j] <- sum(par.y1[j+2,1:K1]*(par.y1[j+2,indyy1]+sum(c(mdm[1:2],unlist(Ds[l,]))* par.y1[j+2,indyy3])))
      
    }
  }
  for(l in 1:nrow(Ds)){
    print(l)
    for(j in 1:(mcmc-2)){
      yyj2[l,j] <- sum(par.y2[j+2,1:K1]*(par.y2[j+2,indyy1]+sum(c(mdm[1:2],unlist(Ds[l,]))* par.y2[j+2,indyy3])))
      
      
    }
  }
  
  
  dgj <- matrix(0,nrow(Ds),(mcmc-2))
  for(i in 1:nrow(Ds)){
    for(j in 1:(mcmc-2)){
      dgj[i,j] <- ifelse(yyj1[i,j]>yyj2[i,j],1,0)
    }
  }
  
  
  return(list(dg21=rowMeans(dg21), dg22=rowMeans(dg22),dgj=rowMeans(dgj)))
}

#Heatmaps for Simulation Setting 1
library(doParallel)  
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores, type="FORK")  
result <- parLapply(cl, sim_1, qlfun_pri)  
stopCluster(cl)  
dgj <- dg21 <- dg22 <- matrix(0,nrow(Ds),length(sim_1))
for(i in 1:length(sim_1)){
  dgj[,i] <- result[[i]]$dgj
  dg21[,i] <- result[[i]]$dg21
  dg22[,i] <- result[[i]]$dg22
}


djt <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dgj))
ggplot(djt, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["1,-1"])) + ylab(expression(C["1,+1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=1/0.8,slope=-1.5/0.8)

datst21 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg21))
ggplot(datst21, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=1/0.8,slope=-1/0.8)

datst22 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg22))
ggplot(datst22, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=1/0.8,slope=-1/0.8)


#Heatmaps for Simulation Setting 2
library(doParallel)  
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores, type="FORK")  
result <- parLapply(cl, sim_2, qlfun_pri)  
stopCluster(cl)  
dgj <- dg21 <- dg22 <- matrix(0,nrow(Ds),length(sim_2))
for(i in 1:length(sim_1)){
  dgj[,i] <- result[[i]]$dgj
  dg21[,i] <- result[[i]]$dg21
  dg22[,i] <- result[[i]]$dg22
}

mdm <- c(-0.5,0,0.5-0.3*0.5)

djt <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dgj))
ggplot(djt, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["1,-1"])) + ylab(expression(C["1,+1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=-0.5*mdm[2]/0.8,slope=1/0.8)

datst21 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg21))
ggplot(datst21, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=0.4*mdm[3]/0.8,slope=1/0.8)

datst22 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg22))
ggplot(datst22, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=0.4*mdm[3]/0.8,slope=1/0.8)


#Heatmaps for Simulation Setting 3
library(doParallel)  
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores, type="FORK")  
result <- parLapply(cl, sim_3, qlfun_pri)  
stopCluster(cl)  
dgj <- dg21 <- dg22 <- matrix(0,nrow(Ds),length(sim_3))
for(i in 1:length(sim_1)){
  dgj[,i] <- result[[i]]$dgj
  dg21[,i] <- result[[i]]$dg21
  dg22[,i] <- result[[i]]$dg22
}

mdm <- c(-0.5,0,0.5-0.3*0.5)

djt <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dgj))
ggplot(djt, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["1,-1"])) + ylab(expression(C["1,+1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=-0.5*mdm[2]/0.8,slope=1/0.8)

datst21 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg21))
ggplot(datst21, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=0.4*mdm[3]/0.8,slope=1/0.8)

datst22 <- data.frame(expand.grid(x=sort(D1su),y=sort(D2su)),z=rowMeans(dg22))
ggplot(datst22, aes(x, y, fill= z)) + geom_tile()+xlab(expression(C["2,+1"])) + ylab(expression(C["2,-1"]))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),legend.position="bottom")+geom_abline(intercept=0.4*mdm[3]/0.8,slope=1/0.8)

