MLECP_MissingData <- function(Nk,N,la){
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
            l1 <- mpfr(10*la,precBits=100)
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
    
    out	
}