#Generation of Simulation Scenarios 1, 2 and 3
#Sample sizes are 500 and 1000

library(copula)
library(truncnorm)
expit <- function(x)exp(x)/(1+exp(x))

#Simulation setting 1
ComplianceSimulate <- function(n, n_sim){
  
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){
    a1 <- 2*rbinom(n,1,.5)-1 #First stage treatment
    
    a2 <- 2*rbinom(n,1,0.5)-1 #Second stage treatment
    
    #Baseline covariates
    X1 <- rnorm(n,-0.5,0.3); u <- rnorm(n,0,1)
    X2 <- rnorm(n,0,sqrt(0.1))
    #Time-varying covariate
    X3<- rnorm(n,0.5+0.3*X1+0.7*X2,0.1)
    
    #Simulate potential compliances from a Gaussian copula with correlation in 'param'
    cop <- normalCopula(param=rep(0.2,6), dim = 4, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=4)
    
    
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm","truncnorm","truncnorm"),
                  paramMargins = list(list(mean=2*u[j],sd=0.5,a=0,b=1),
                                      list(mean=1-u[j],sd=0.5,a=0,b=1),
                                      list(mean=0.5*u[j],sd=0.5,a=0,b=1),
                                      list(mean=1.5*u[j],sd=0.5,a=0,b=1)))
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D1 <- D[,1] #Stage 1: A1=+1
    D2 <- D[,2] #Stage 1: A1=-1 
    D3 <- D[,3] #Stage 2: A2 =+1
    D4 <- D[,4]
    s <- rep(NA,n)
    
    #First stage response indicators
    s[a1==1] <- rbinom(length(which(a1==1)),1,exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))/(1+exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))))
    s[a1==-1] <- rbinom(length(which(a1==-1)),1,exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))/(1+exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))))
    
    #############
    #Simulate outcomes without noise
    
    y1 <- rep(NA,n)
    
   y1[a1==1&s==1] = 0.7+0.2*X1[a1==1&s==1]+0.2*(X2[a1==1&s==1]+D2[a1==1&s==1])+
     2*(-1+1.5*D1[a1==1&s==1]+0.8*D2[a1==1&s==1])
    
    y1[a1==1&s==0&a2==1] = 0.7+0.2*X1[a1==1&s==0&a2==1]+0.2*X2[a1==1&s==0&a2==1]+
      0.2*D2[a1==1&s==0&a2==1]+2*(-1+1.5*D1[a1==1&s==0&a2==1]+0.8*D2[a1==1&s==0&a2==1])+
      0.5*X3[a1==1&s==0&a2==1]+0.5*D4[a1==1&s==0&a2==1]+2*(-1+D3[a1==1&s==0&a2==1]+0.8*D4[a1==1&s==0&a2==1])
    
    y1[a1==1&s==0&a2==-1] = 0.7+0.2*X1[a1==1&s==0&a2==-1]+0.2*X2[a1==1&s==0&a2==-1]+
      0.2*D2[a1==1&s==0&a2==-1]+2*(-1+1.5*D1[a1==1&s==0&a2==-1]+0.8*D2[a1==1&s==0&a2==-1])+
      0.5*X3[a1==1&s==0&a2==-1]+0.5*D4[a1==1&s==0&a2==-1]-2*(-1+D3[a1==1&s==0&a2==-1]+0.8*D4[a1==1&s==0&a2==-1])
    
    y1[a1==-1&s==1] = 0.7+0.2*X1[a1==-1&s==1]+0.2*(X2[a1==-1&s==1]+D2[a1==-1&s==1])-2*(-1+1.5*D1[a1==-1&s==1]+0.8*D2[a1==-1&s==1])
    
    y1[a1==-1&s==0&a2==1] = 0.7+0.2*X1[a1==-1&s==0&a2==1]+0.2*(X2[a1==-1&s==0&a2==1]+D2[a1==-1&s==0&a2==1])-2*(-1+1.5*D1[a1==-1&s==0&a2==1]+0.8*D2[a1==-1&s==0&a2==1])+
      0.5*X3[a1==-1&s==0&a2==1]+0.5*D4[a1==-1&s==0&a2==1]+
      2*(-1+D3[a1==-1&s==0&a2==1]+0.8*D4[a1==-1&s==0&a2==1])
    
    y1[a1==-1&s==0&a2==-1] = 0.7+0.2*X1[a1==-1&s==0&a2==-1]+0.2*(X2[a1==-1&s==0&a2==-1]+
                                                                   D2[a1==-1&s==0&a2==-1])-2*(-1+1.5*D1[a1==-1&s==0&a2==-1]+0.8*D2[a1==-1&s==0&a2==-1])+
      +0.5*X3[a1==-1&s==0&a2==-1]+0.5*D4[a1==-1&s==0&a2==-1]-2*(-1+D3[a1==-1&s==0&a2==-1]+0.8*D4[a1==-1&s==0&a2==-1])
    
    
    #Add noise to form outcomes on the scale of the coefficients
    y <- y1 + rnorm(n,0,0.1) 
    
    sim_list[[i]] <- data.frame(a1,a2,D1,D2,D3,D4,s,y,y1,X1,X2,X3,u)
    
    print(i)
  }
  sim_list
}

#set.seed(38265)
sim_1 <- ComplianceSimulate(n=1000,n_sim=1)


#Simulation setting 2
ComplianceSimulate <- function(n, n_sim){
  
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){
    a1 <- 2*rbinom(n,1,.5)-1 #First stage treatment
    
    a2 <- 2*rbinom(n,1,0.5)-1 #Second stage treatment
    
    #Baseline covariates
    X1 <- rnorm(n,-0.5,0.3); u <- rnorm(n,0,1)
    X2 <- rnorm(n,0,sqrt(0.1))
    #Time-varying covariate
    X3<- rnorm(n,0.5+0.3*X1+0.7*X2,0.1)
    
    #Simulate potential compliances from a Gaussian copula with correlation in 'param'
    cop <- normalCopula(param=rep(0.2,6), dim = 4, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=4)
    
    
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm","truncnorm","truncnorm"),
                  paramMargins = list(list(mean=2*u[j],sd=0.5,a=0,b=1),
                                      list(mean=1-u[j],sd=0.5,a=0,b=1),
                                      list(mean=0.5*u[j],sd=0.5,a=0,b=1),
                                      list(mean=1.5*u[j],sd=0.5,a=0,b=1)))
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D1 <- D[,1] #Stage 1: A1=+1
    D2 <- D[,2] #Stage 1: A1=-1 
    D3 <- D[,3] #Stage 2: A2 =+1
    D4 <- D[,4]
    s <- rep(NA,n)
    
    #First stage response indicators
    s[a1==1] <- rbinom(length(which(a1==1)),1,exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))/(1+exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))))
    s[a1==-1] <- rbinom(length(which(a1==-1)),1,exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))/(1+exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))))
    
    #############
    #Simulate outcomes without noise
    
    y1 <- rep(NA,n)
    
    y1[a1==1&s==1] = 0.7+0.2*X1[a1==1&s==1]+0.2*(X2[a1==1&s==1]+D2[a1==1&s==1])+
      2*(1.5*D1[a1==1&s==1]-0.8*D2[a1==1&s==1]-0.5*X2[a1==1&s==1])
    
    y1[a1==1&s==0&a2==1] = 0.7+0.2*X1[a1==1&s==0&a2==1]+0.2*X2[a1==1&s==0&a2==1]+
      0.2*D2[a1==1&s==0&a2==1]+2*(1.5*D1[a1==1&s==0&a2==1]-0.8*D2[a1==1&s==0&a2==1]-0.5*X2[a1==1&s==0&a2==1])+
      0.5*X3[a1==1&s==0&a2==1]+0.5*D4[a1==1&s==0&a2==1]+2*(D3[a1==1&s==0&a2==1]-0.8*D4[a1==1&s==0&a2==1]+0.4*X3[a1==1&s==0&a2==1])
    
    y1[a1==1&s==0&a2==-1] = 0.7+0.2*X1[a1==1&s==0&a2==-1]+0.2*X2[a1==1&s==0&a2==-1]+
      0.2*D2[a1==1&s==0&a2==-1]+2*(1.5*D1[a1==1&s==0&a2==-1]-0.8*D2[a1==1&s==0&a2==-1]-0.5*X2[a1==1&s==0&a2==-1])+
      0.5*X3[a1==1&s==0&a2==-1]+0.5*D4[a1==1&s==0&a2==-1]-2*(D3[a1==1&s==0&a2==-1]-0.8*D4[a1==1&s==0&a2==-1]+0.4*X3[a1==1&s==0&a2==-1])
    
    y1[a1==-1&s==1] = 0.7+0.2*X1[a1==-1&s==1]+0.2*(X2[a1==-1&s==1]+D2[a1==-1&s==1])-2*(1.5*D1[a1==-1&s==1]-0.8*D2[a1==-1&s==1]-0.5*X2[a1==-1&s==1])
    
    y1[a1==-1&s==0&a2==1] = 0.7+0.2*X1[a1==-1&s==0&a2==1]+0.2*(X2[a1==-1&s==0&a2==1]+D2[a1==-1&s==0&a2==1])-
      2*(1.5*D1[a1==-1&s==0&a2==1]-0.8*D2[a1==-1&s==0&a2==1]-0.5*X2[a1==-1&s==0&a2==1])+
      0.5*X3[a1==-1&s==0&a2==1]+0.5*D4[a1==-1&s==0&a2==1]+
      2*(D3[a1==-1&s==0&a2==1]-0.8*D4[a1==-1&s==0&a2==1]+0.4*X3[a1==-1&s==0&a2==1])
    
    y1[a1==-1&s==0&a2==-1] = 0.7+0.2*X1[a1==-1&s==0&a2==-1]+0.2*(X2[a1==-1&s==0&a2==-1]+D2[a1==-1&s==0&a2==-1])-
      2*(1.5*D1[a1==-1&s==0&a2==-1]-0.8*D2[a1==-1&s==0&a2==-1]-0.5*X2[a1==-1&s==0&a2==-1])+
      +0.5*X3[a1==-1&s==0&a2==-1]+0.5*D4[a1==-1&s==0&a2==-1]-2*(D3[a1==-1&s==0&a2==-1]-0.8*D4[a1==-1&s==0&a2==-1]+0.4*X3[a1==-1&s==0&a2==-1])
    
    
    #Add noise to form outcomes on the scale of the coefficients
    y <- y1 + rnorm(n,0,0.1) 
    
    sim_list[[i]] <- data.frame(a1,a2,D1,D2,D3,D4,s,y,y1,X1,X2,X3,u)
    
    print(i)
  }
  sim_list
}

#set.seed(38265)
sim_2 <- ComplianceSimulate(n=1000,n_sim=1)

#Simulation setting 3
omplianceSimulate <- function(n, n_sim){
  
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){
    a1 <- 2*rbinom(n,1,.5)-1 #First stage treatment
    
    a2 <- 2*rbinom(n,1,0.5)-1 #Second stage treatment
    
    #Baseline covariates
    X1 <- rnorm(n,-0.5,0.3); u <- rnorm(n,0,1)
    X2 <- rnorm(n,0,sqrt(0.1))
    #Time-varying covariate
    X3<- rnorm(n,0.5+0.3*X1+0.7*X2,0.1)
    
    #Simulate potential compliances from a Gaussian copula with correlation in 'param'
    cop <- normalCopula(param=rep(0.2,6), dim = 4, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=4)
    
    
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm","truncnorm","truncnorm"),
                  paramMargins = list(list(mean=2*u[j]-0.2*X1[j],sd=0.5,a=0,b=1),
                                      list(mean=1-u[j],sd=0.5,a=0,b=1),
                                      list(mean=0.5*u[j]+0.5*X2[j],sd=0.5,a=0,b=1),
                                      list(mean=1.5*u[j]+0.3*X3[j],sd=0.5,a=0,b=1)))
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D1 <- D[,1] #Stage 1: A1=+1
    D2 <- D[,2] #Stage 1: A1=-1 
    D3 <- D[,3] #Stage 2: A2 =+1
    D4 <- D[,4]
    s <- rep(NA,n)
    
    #First stage response indicators
    s[a1==1] <- rbinom(length(which(a1==1)),1,exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))/(1+exp(1*D1[a1==1]-1.5+0.2*(X2[which(a1==1)]))))
    s[a1==-1] <- rbinom(length(which(a1==-1)),1,exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))/(1+exp(1*D2[a1==-1]-1.5+0.3*(X1[which(a1==-1)]))))
    
    #############
    #Simulate outcomes without noise
    
    y1 <- rep(NA,n)
    
    y1[a1==1&s==1] = 0.7+0.2*X1[a1==1&s==1]+0.2*(X2[a1==1&s==1]+D2[a1==1&s==1])+
      2*(1.5*D1[a1==1&s==1]-0.8*D2[a1==1&s==1]-0.5*X2[a1==1&s==1])
    
    y1[a1==1&s==0&a2==1] = 0.7+0.2*X1[a1==1&s==0&a2==1]+0.2*X2[a1==1&s==0&a2==1]+
      0.2*D2[a1==1&s==0&a2==1]+2*(1.5*D1[a1==1&s==0&a2==1]-0.8*D2[a1==1&s==0&a2==1]-0.5*X2[a1==1&s==0&a2==1])+
      0.5*X3[a1==1&s==0&a2==1]+0.5*D4[a1==1&s==0&a2==1]+2*(D3[a1==1&s==0&a2==1]-0.8*D4[a1==1&s==0&a2==1]+0.4*X3[a1==1&s==0&a2==1])
    
    y1[a1==1&s==0&a2==-1] = 0.7+0.2*X1[a1==1&s==0&a2==-1]+0.2*X2[a1==1&s==0&a2==-1]+
      0.2*D2[a1==1&s==0&a2==-1]+2*(1.5*D1[a1==1&s==0&a2==-1]-0.8*D2[a1==1&s==0&a2==-1]-0.5*X2[a1==1&s==0&a2==-1])+
      0.5*X3[a1==1&s==0&a2==-1]+0.5*D4[a1==1&s==0&a2==-1]-2*(D3[a1==1&s==0&a2==-1]-0.8*D4[a1==1&s==0&a2==-1]+0.4*X3[a1==1&s==0&a2==-1])
    
    y1[a1==-1&s==1] = 0.7+0.2*X1[a1==-1&s==1]+0.2*(X2[a1==-1&s==1]+D2[a1==-1&s==1])-2*(1.5*D1[a1==-1&s==1]-0.8*D2[a1==-1&s==1]-0.5*X2[a1==-1&s==1])
    
    y1[a1==-1&s==0&a2==1] = 0.7+0.2*X1[a1==-1&s==0&a2==1]+0.2*(X2[a1==-1&s==0&a2==1]+D2[a1==-1&s==0&a2==1])-
      2*(1.5*D1[a1==-1&s==0&a2==1]-0.8*D2[a1==-1&s==0&a2==1]-0.5*X2[a1==-1&s==0&a2==1])+
      0.5*X3[a1==-1&s==0&a2==1]+0.5*D4[a1==-1&s==0&a2==1]+
      2*(D3[a1==-1&s==0&a2==1]-0.8*D4[a1==-1&s==0&a2==1]+0.4*X3[a1==-1&s==0&a2==1])
    
    y1[a1==-1&s==0&a2==-1] = 0.7+0.2*X1[a1==-1&s==0&a2==-1]+0.2*(X2[a1==-1&s==0&a2==-1]+D2[a1==-1&s==0&a2==-1])-
      2*(1.5*D1[a1==-1&s==0&a2==-1]-0.8*D2[a1==-1&s==0&a2==-1]-0.5*X2[a1==-1&s==0&a2==-1])+
      +0.5*X3[a1==-1&s==0&a2==-1]+0.5*D4[a1==-1&s==0&a2==-1]-2*(D3[a1==-1&s==0&a2==-1]-0.8*D4[a1==-1&s==0&a2==-1]+0.4*X3[a1==-1&s==0&a2==-1])
    
    
    #Add noise to form outcomes on the scale of the coefficients
    y <- y1 + rnorm(n,0,0.1) 
    
    sim_list[[i]] <- data.frame(a1,a2,D1,D2,D3,D4,s,y,y1,X1,X2,X3,u)
    
    print(i)
  }
  sim_list
}
sim_3 <- ComplianceSimulate(n=1000,n_sim=1)
