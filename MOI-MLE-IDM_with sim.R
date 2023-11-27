# Title        : Script for the MLE of lineage frequencies and MOI parameter using Incomplete Data Model
# Objective    : Estimate lineage frequencies, MOI,
# Created by   : Meraj Hashemi, Kristan A. Schneider
# Created on   : 25.04.22
# Last modified: 07.03.23


library('ggplot2')
library('stringr')

################################# This function imports data #################################

DatImp<-function(path){
  if(substring(path,nchar(path)-3,nchar(path))==".xls"){
    dat <- openxlsx::read.xlsx(path,1)
  }
  else{
    if(substring(path,nchar(path)-4,nchar(path))==".xlsx"){
      dat <- openxlsx::read.xlsx(path,1)
    }
    else{
      if(substring(path,nchar(path)-3,nchar(path))==".txt"){
        dat <- read.table(path,header=TRUE, sep="\t")
      }
      else{
        if(substring(path,nchar(path)-3,nchar(path))==".csv"){
          dat <- read.csv(path,header=TRUE,sep=";")
        }
      }
    }  
  }
  dat
}     

dat <- DatImp("/Users/kristanschneider/Library/CloudStorage/OneDrive-Persönlich/Documents/Forschung/Likelihood robustness/Submission/plos-latex-template/Final/Example Data/STR.xlsx")
iii <-Nk(dat)

MLE_OM(iii[[1]],iii[[2]],iii[[3]])
################################ This function calculate Nk #################################

Nk <- function(dat){
  for(k in 1:nrow(dat)){
    if(dat[k,1]==""||is.na(dat[k,1])){
      dat[k,1] <- dat[k-1,1]
    }
  }
  N <- length(unique(dat[,1]))
  dat <- dat[!is.na(dat[,2]),]
  Nplus <- length(unique(dat[,1]))
  Nknum <- length(unique(dat[,2]))
  dat <- dat[!duplicated(dat),]
  out <- list(N,t(as.matrix(table(dat[,2]))),N-Nplus)
  names(out) <- c("N","N_k","n_0")
  out
}

########################################################################################
#------------------------- Functions for calculating the MLEs---------------------------
########################################################################################

MLE <- function(N, N_k, n_0=0, model = "IDM", lambda_initial = 1, eps_initial=0.1){
  if(model == "OM"){
    print(n_0)
    MLE_OM(N-n_0,N_k,lambda_initial)
  }else if(model == "IDM"){
      MLE_IDM(N, N_k,n_0,lambda_initial, eps_initial)
  }else{
    Warning("option model needs to be eiter `IDM' or `OM'")
  }
}    

################################ The MLE of the IDM #################################

#' The EM algorithm to derive the MLE for the IDM
#'
#' @param N integer; Sample size
#' @param N0 integer; number of empty records
#' @param Nk vector of integers;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param lambda_initial float; initial value of lambda
#' @param eps_initial float; inital value of epsilon (probability of the
#'   lineages remain undetected)
#'
#' @return the function returns a list of values as follows:
#'         1) the MLE of the probability of lineages remaining undetected (epsilon)
#'         2) the MLE of the MOI parameter (lambda)
#'         3) the MLE of the average MOI (psi)
#'         4) the MLE of the lineage frequencies
#'
#' @examples EM(40, 1, c(23,27), 1, 0.1)
MLE_IDM <- function(N, Nk,n0, lambda_initial, eps_initial) {
  thr1 <- 10^-8
  thr2 <- 10^-12
  thr3 <- 10^-20
  z <- Nk
  Nk <- Nk[Nk!=0]
  n <- length(Nk)
  Nnk <- N - Nk
  snk <- sum(Nk)
  #initial values
  lamt <- 0.1
  lamnext <- lambda_initial
  pkt <- as.vector(array(1/n,c(1,n))) + 0.1
  pnext <- as.vector(array(1/n,c(1,n)))
  epst <- 0.01
  epsnext <- eps_initial 
  while(abs(lamt - lamnext) +abs(epst - epsnext) + sqrt(sum((pkt - pnext)^2)) > thr1) {
    lamt <- lamnext
    epst <- epsnext
    pkt <- pnext
    nextiter <- EM_next_iteration(lamnext, pnext, epst, N, n0, Nk, Nnk, snk)
    pnext <- nextiter[[1]]
    epsnext<- nextiter[[2]] 
    wt <- nextiter[[3]]
    ntt <- wt/N
    lamnext <-  ntt + 1
    lamnext <- EM_lambda_Newton(lamnext, N, thr1, ntt)
  }
  
  pnextz <- array(0,length(z))  
  pnextz[z > 0] <- pnext 
  psi <- lamnext/(1-exp(-lamnext))
  final <- list(epsnext, lamnext, psi, pnextz)
  #names(final) <- c("")
  final
}

#---------------------------------internal function-------------------------------

EM_lambda_Newton <- function(initial, N, thr, ntt) {
  lamt <- 0
  lamnext <- initial
  while (abs(lamt - lamnext) > thr) {
    lamt <- lamnext
    exp_l <- 1 - exp(-lamt)
    newt <- 1 - exp(-lamt)*ntt
    lamnext <- (ntt*(1 - exp(-lamt)*(1 + lamt)))/newt 
  }
  
  lamnext
}

#---------------------------------internal function-------------------------------

EM_next_iteration <- function (lamnext, pkt, epst, N, N0, Nk, Nnk, snk) {
  #prereuisite
  expt <- exp(lamnext*pkt)
  exp1t <- expt - 1
  exp1expt <- expt/exp1t
  expepst <- epst*exp1t + 1
  if (N0 > 0) {
    tt <- N0/(-1 + prod(expepst))
  }
  else {
    tt <- 0
  }
  
  exp2epst <- expt/expepst
  
  wt <- lamnext*( sum(pkt*(Nk*exp1expt + Nnk*epst*exp2epst)) 
                  + epst*tt*sum(pkt*exp2epst))
  exp1e <- exp1t/expepst
  vkt <- Nnk*exp1e
  vt <- epst*(sum(vkt) + sum(exp1e)*tt)
  ukt <- pkt*lamnext*(Nk*exp1expt + epst*exp2epst*(Nnk + tt))
  #next iteration
  pnext <- ukt/sum(ukt)
  epsnext <- 1/(1 + (snk/vt))
  list(pnext, epsnext, wt)
}

################################ The MLE of the OM #################################

#' function to derive the MLE for the original model
#'
#' @param N integer; Sample size
#' @param Nk vector of integers;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param la float; initial value of lambda
#'
#' @return the function returns a list of values as follows: 1) the MLE of the
#'   MOI parameter (lambda) 2) the MLE of the average MOI (psi) 3) the MLE of
#'   the lineage frequencies
#' @export
#'
#' @examples MLE(c(22,25,49,32,18),97)
MLE_OM <- function(N,Nk,la=1){
  sel <- Nk
  Nk <- sel[sel>0]
  nk <- Nk/N
  l1 <- 2.5         # initial value
  l0 <- 0
  eps <- 10^(-8)       # precision 
  k <- 1
  while(abs(l0-l1)>eps && k<50 && l1>0){
    k <- k+1
    l0 <- l1
    l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
  }
  if(k==50 || l1<0){
    for(st in 1:10){
      l1 <- st
      l0 <- l1+1
      k <- 1
      while(abs(l0-l1)>eps && k<100 && l1>0){
        k <- k+1
        l0 <- l1
        l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
      }
      if(abs(l0-l1)<eps){
        break
      }
    }
    if(abs(l0-l1)>eps){
      l1 <- Rmpfr::mpfr(10*la,precBits=100)
      l0 <- l1+1
      while(abs(l0-l1)>eps){
        l0 <- l1
        l1=l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
      }
    }        
  }
  pk <- -1/l1*log(1-nk*(1-exp(-l1))) 
  pk1 <- array(0,length(sel))  
  pk1[sel>0] <- pk  
  psi <- l1/(1-exp(-l1))
  out <- list(l1,psi,pk1)
  names(out) <- c("MOI parameter lambda","average MOI","lineage frequencies")
  out	
}

MLE(97,c(22,25,49,32,18),model="OM")

########################################################################################
#------------------------- Functions for the simulation study --------------------------
########################################################################################

####################################### cpoiss #########################################

#' Generates conditional Poisson random numbers
#'
#' @param lambda float; the MOI parameter
#' @param N integer; the sample size 
#'
#' @return a vector of randomly generated conditional Poisson numbers
#'
#' @examples  cpoiss(1.5, 10)
#' 
cpoiss<-function(lambda,N){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- rep(0,N)
  x <- runif(N,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:N){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k==m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out[i] <- k
  }
  out
}

####################################### mnom #########################################

#' Generates molecular dataset for a given set of model parameters
#'
#' @param M either a positive integer or a vector of positive integers
#'   corresponding to conditional Poisson random numbers
#' @param p vector; vector of lineage frequencies
#'
#' @return a 0-1 matrix of size N x n where each row corresponds to a sample
#'   with
#'
#' @examples 
#'  mnom(8, c(0.25,0.25,0.25,0.25))
#'  
#'  mnom(c(8,5,6), c(0.25,0.25,0.25,0.25))
#' 
mnom <- function(M, p){
  N <- length(M)
  out<-matrix(0, N, length(p))
  for(k in 1:N){
    out[k,] <- rmultinom(1,M[k],p)
  }
  out <- out
  out
}

####################################### runsim #########################################

#' Generates a molecular dataset with incomplete information 
#'
#' @param data matrix; a 0-1 matrix corresponding to N blood samples 
#' @param eps float; the probability of lineages remaining undetected
#' @param N integer; sample size
#' @param n integer; number of lineages
#'
#' @return a list with the following values:
#'            1) number of empty records
#'            2) the vector of observed prevalences
#'            3) dataset with incomplete information
#'
#' @examples
IncompleteData <- function(data, eps, N, n){
  ran <- runif(N*n)
  ran <- matrix((ran > eps)*1, N, n)
  ndata <- ran*data
  nkmiss <- colSums(ndata)
  nzero <- sum((rowSums(ndata) == 0)*1)
  list(nzero,nkmiss,data)
}


####################################### runsim #########################################


#' Runs a simulation
#'
#' @param K number of simulation steps
#' @param eps vector of values for probability of lineages remaining undetected,
#'   e.g., eps <- c(0.05, 0.1, 0.15)
#' @param lambda vector of MOI parameter values, e.g., 
#' lambda <- seq(0.1, 2.5, 0.1)
#' @param linfreq list; list of lineage-frequency distributions, e.g., 
#' linfreq <- list(c(0.5, 0.5), c(0.75, 0.25))
#' @param ssize vector of sample size values, e.g., 
#' ssize <- sort(c(50,100,150,200), decreasing = TRUE)
#' @param path the path to store the simulation results
#'
#' @return a .txt-file containing the simulation results for each
#'   lineage-frequency distribution. The file containing the simulation results
#'   corresponding to a lineage-frequency distribution p with n lineages is
#'   stored under the given path in a .txt-file which is named as follows:
#'   data-n<n>-maxfreq<max(p)>.txt. A line of the file contains the following
#'   results which are stored with an order exactly as mentioned below: 
#'   1) true MOI parameter (lambda) 
#'   2) true probability of lineages remain undetected(varepsilon) 
#'   3) the empirical mean of the following variables based on the
#'   K-simulated datasets:
#'              - the MLE of varepsilon
#'              - the MLE of lambda from IDM
#'              - the MLE of lambda from OM
#'              - the MLE of average MOI, i.e., psi from IDM
#'              - the MLE of average MOI, from OM
#'              - the MLE of lineage frequencies from IDM
#'              - the MLE of lineage frequencies from OM
#'              - the Euclidean norm of frequency distribution from IDM
#'              - the Euclidean norm of frequency distribution from OM
#'              - true prevalence from the IDM
#'              - observed prevalence from the IDM
#'              - true prevalence from the OM
#'   4) the empirical variance of the above variables based on the K-simulated
#'   datasets 
#'   5) sample size 
#'   6) lineage-frequency distribution
#' @examples
runsim <- function(K, eps, lambda, linfreq, ssize, path){
  s <- Sys.time()
  ssize <- sort(ssize, decreasing = TRUE)
  NN <- ssize[1]
  restN <- ssize[-1]
  sz <- length(ssize)
  for (p in linfreq) {
    pp <- 100*round(p[1],digits = 2)
    n <- length(p)
    for (e in eps) {
      for (lam in lambda) {
        simfinal <- rep(list(array(NA, c(K, 7 + 5*n))), sz)
        sim <- 0
        dgen <- rep(list(list(0)), sz)
        while (sim < K) {
          sim <- sim + 1
          md <- innersamplegenerator_MissingData(NA, NA, NN, NA, e, lam, p, n)
          N0 <- md[[1]]
          Nk <- md[[2]]
          data <- md[[3]]
          em <- EM(NN, N0, Nk, lam, e)
          mle <- MLE(Nk, NN - N0, lam)
          
          enem <- sqrt(sum((em[[4]]-p)^2))  ## Euclidian norm freqs - true parameter  - extended model
          enom <- sqrt(sum((mle[[3]]-p)^2))  ## Euclidian norm freqs - true parameter - original model
          
          elpk <- exp(-em[[2]] * em[[4]])
          tpem <- (1-elpk)/(1-exp(-em[[2]]))                              ## true prevalence  - extended model
          opem <- (1-em[[1]])*(1-elpk)/(1-prod(em[[1]]*(1-elpk)+elpk))    ## observed prevalence  -   extended model
          tpom <- (1-exp(-mle[[1]] * mle[[3]]))/(1-exp(-mle[[1]]))        ## true prevalence  - original model
          
          # ^eps, ^lam, ^lam*, ^psi, ^psi*, ^pp, ^pp* 
          result <- c(em[[1]], em[[2]], mle[[1]], em[[3]], mle[[2]], em[[4]], mle[[3]],enem,enom,tpem,opem,tpom)
          
          
          simfinal[[1]][sim,] <- result
          t <- 1
          for (N in restN) {
            #print(N)
            t <- t + 1
            dataN <- data[1:N,]
            mdN <- IncompleteData(dataN, e, N, n)
            N0N <- mdN[[1]]
            NkN <- mdN[[2]]
            mdN <- innersamplegenerator_MissingData(N0N, NkN, N, dataN, e, lam, p, n)
            N0N <- mdN[[1]]
            NkN <- mdN[[2]]
            data <- md[[3]]
            em <- EM(N, N0N, NkN, lam, e)
            mle <- MLE(NkN, N - N0N, lam)
            
            enem <- sqrt(sum((em[[4]]-p)^2))  ## Euclidian norm freqs - true parameter  - extended model
            enom <- sqrt(sum((mle[[3]]-p)^2))  ## Euclidian norm freqs - true parameter - original model
            
            elpk <- exp(-em[[2]] * em[[4]])
            tpem <- (1-elpk)/(1-exp(-em[[2]]))                              ## true prevalence  - extended model
            opem <- (1-em[[1]])*(1-elpk)/(1-prod(em[[1]]*(1-elpk)+elpk))    ## observed prevalence  -   extended model
            tpom <- (1-exp(-mle[[1]] * mle[[3]]))/(1-exp(-mle[[1]]))        ## true prevalence  - original model
            
            # ^eps, ^lam, ^lam*, ^psi, ^psi*, ^pp, ^pp* 
            result <- c(em[[1]], em[[2]], mle[[1]], em[[3]], mle[[2]], em[[4]], mle[[3]],enem,enom,tpem,opem,tpom)
            
            simfinal[[t]][sim,] <-  result
          }
        }
        
        for(d in 1:sz){
          final <- c(lam, e, apply(simfinal[[d]],2,mean), apply(simfinal[[d]],2,var), ssize[d], p)
          # writting the results to the file
          write.table(t(final),paste(path,"/data-n",toString(n), "-maxfreq",pp, ".txt",sep=""),
                      append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
        }
      }
    }
  }
  print(Sys.time() - s)
}


#---------------------------------internal function-------------------------------
#
###for nested sample generator, generates sample of size N
innersamplegenerator_MissingData <- function(N0,Nk,N,data,e,lambda,p,n) {
  regular <- FALSE
  Nn <- 1/N^n
  if (is.na(N0) == TRUE) {
    M <- cpoiss(lambda, N)
    data <- sign(mnom(M, p))
    result <- IncompleteData(data, e, N, n)
    N0 <- result[[1]]
    Nk <- result[[2]]
    data <- result[[3]]
  }
  if (sum(Nk > 0) == 1) {
    regular <- FALSE
  }
  else if (N0 != 0 & (sum(Nk) == (N - N0) || prod(1 - Nk/N) - N0/N < Nn )){  ### prod(1 - Nk/N) - N0/N == 0 can lead to numerical errors e.g., prod(1-c(30,30)/50)==8/50 gives FALSE, difference 2.775558e-17.   
    ## hence abs(prod(1 - Nk/N) - N0/N) < 1/(N^n) gets artifacts from numerical precision
    regular <- FALSE
  }
  else if (N0 == 0 & (sum(Nk) == N || max(Nk) == N)){
    regular <- FALSE
  }
  else {
    regular <- TRUE
  }
  while (regular == FALSE){
    M <- cpoiss(lambda, N)
    data <- sign(mnom(M, p))
    result <- IncompleteData(data, e, N, n)
    N0 <- result[[1]]
    Nk <- result[[2]]
    if (sum(Nk > 0) == 1) {
      regular <- FALSE
    }
    else if (N0 != 0 & (sum(Nk) == (N - N0) || prod(1 - Nk/N) - N0/N <Nn)){
      regular <- FALSE
    }
    else if (N0 == 0 & (sum(Nk) == N || max(Nk) == N)){
      regular <- FALSE
    }
    else {
      regular <- TRUE
    }
  } 
  list(N0,Nk,data)
}

