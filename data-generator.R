
########################################################################################
#--------------------------Data Generator---------------------------------------------
########################################################################################

#Generates conditional Poisson random numbers
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
#------------------------------------------------------------------------------------
###Generates molecular dataset for model parameters
mnom <- function(N,lambda,p,n){
   
    cp<-cpoiss(lambda,N)
    out<-matrix(0,N,n)
    for(k in 1:N){
        out[k,] <- rmultinom(1,cp[k],p)
    }
    out <- sign(out)
    out
}
#------------------------------------------------------------------------------------
###Generates multinomial random numbers based on conditional Poisson random numbers 
zero_one_matrix <- function (cp, N, n, p) {
    out<-matrix(0,N,n)
    for(k in 1:N){
        out[k,] <- rmultinom(1,cp[k],p)
    }
    out <- sign(out)
    out
}





########################################################################################
#-----------------------------Negative Binomial---------------------------------------
########################################################################################

###generates nbinom numbers
### 0<mu<100
c_negb<-function(N, success, p){
    m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
    out <- rep(0,N)
    x <- runif(N,min=0,max=1)
    p0 <- pnbinom(0,size = success, prob =  p)
    nc <- 1/(1 - p0)
    pvec <- (pnbinom(1:m,size = success, prob = p) - p0)*nc
    pvec <- c(pvec,1)
    for (i in 1:N){
        k <- 1
        while(x[i] > pvec[k]){
            k <- k+1
        }
        if(k==m){ # if a m>=100 is drawn this is executed
            k <- k+1
            a <- dnbinom(k, size = success, prob = p)*nc
            b <- pvec[m]+a
            while(x[i]>b){
                k <- k+1
                a <- a*(1-p)*(k+success-1)/k
                b <- b+a
            }
        }
        out[i] <- k
    }
    out
}
#------------------------------------------------------------------------------------


###generates random numbers based on nbinom model
m_negb<-function(N,success,p,pk,n){
    nb<-c_negb(N, success, p)
    out<-matrix(0,N,n)
    for(k in 1:N){
        out[k,]=rmultinom(1,nb[k],pk)
    }
    sign(out)
}


#------------------------------------------------------------------------------------



###findes values corresponding to different lambda and dispersion factors
equ_lam <- function(lambda,alpha){
    k <- 0
    out <- list()
    for(lam in lambda){
        for(a in alpha){
            si <- lam/(1-exp(-lam))
            extr <-  si + a*(lam + 1 - si)
            thr <- 10^(-3)
            p0 <- 1
            p <- 0.9
            k <- k + 1
            while(abs(p0 - p)>thr){
                p0 <- p
                fp <- (p*extr - 1)/(1-p)
                hp <-  (extr - 1/p)*(1/(1-p^fp)) - si
                dhp <- p^(-2)/(1 - p^fp) + (extr - 1/p)*(p^fp)*((extr - 1)*log(p)/((1-p)^2) + fp/p)/((1-p^fp)^2)
                p <- p - hp/dhp
            }
            
            success <- (p*extr - 1)/(1-p)
            out[[k]] <- c(lam, a, success, p)
        }
    }
    out <- t(matrix(unlist(out), 4, length(out)))
    out
}
#------------------------------------------------------------------------------------




