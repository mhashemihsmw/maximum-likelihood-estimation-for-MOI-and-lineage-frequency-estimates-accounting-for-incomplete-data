p <- 1:2
p <- p/sum(p)
p <- c(0.75,0.15,0.05,0.05)
#p <- c(0.5,0.5)
p <- c(1,1,1,1)/4
eps <- 0.15
lam <- 2
N <- 50
out <-c()


FI <- function (N,lam,p,eps=0) {
  #prerequisite
  n <-length(p)
  if(eps==0){
    FI <- array(0,c(n+2,n+2))
    iFI <- array(0,c(n+2,n+2))
    el <- exp(lam)
    eml <-1/el
    elk <- exp(lam*p)
    eplkmo <- (elk-1)
    
    
    ## d^2L/dl^2
    FI[1,1] <- -1/(1-eml)+el*sum( p^2/eplkmo)
    
    ## d^2L/dldpk
    
    FI[1,1:n+2] <- - el*(1- p/eplkmo*lam )
    FI[1:n+2,1] <- FI[1,1:n+2]
    
    ## d^2L/dpk^2
    
    D <- (el*lam^2/eplkmo)
    #FI <- FI + diag(c(0,0,D ))
    
    FI <- N/(el-1)*FI
    FI[1:n+2,2] <- 1
    FI[2,1:n+2] <- FI[1:n+2,2] 
    
    ## use blockwise inverision
    D <- N/(el-1)*D
    D1 <- 1/D
    d <- -sum(D1)
    D2 <- c(-1,D1)
    D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
    
    A0 <-FI[1,1]
    B <- FI[1,2:(n+2)]
    A1 <- FI[1:n+2,1:n+2]
    DC <- D3 %*% B
    A1 <- solve(A0 - (B %*% DC))
    
    iFI[1,1] <- A1
    C2 <-  -DC %*% A1 
    iFI[0:n+2,1] <- C2
    B2 <- t(C2)
    iFI[1,0:n+2] <- B2
    iFI[0:n+2,0:n+2] <- D3 - DC %*% B2 
    iFI <- iFI[-2,-2]
    
  }else{
    FI <- array(0,c(n+3,n+3))
    iFI <- array(0,c(n+3,n+3))
    el <- exp(lam)
    eml <-1/el
    elk <- exp(lam*p)
    Ak <- (elk-1)*eps +1
    tau <- prod(Ak)
    tau <- tau/(tau-1)
    pkAk <- t(p/Ak)
    eplkmo <- (elk-1)
    eplmoAl <- t(eplkmo/Ak)
    Tl <- eps * p * elk /Ak
    Te <- eplkmo/Ak
    Tp <- eps*lam*elk/Ak
    ## d^2L/dl^2
    FI[1,1] <- -1/(1-eml)+(1-eps)*el*(pkAk%*% (p/eplkmo)) + tau*sum(Tl)^2
    
    ## d^2L/dlde
    FI[1,2] <- -el*sum(pkAk)  + tau * sum(Tl)*sum(Te)
    FI[2,1] <- FI[1,2]
    
    ## d^2L/dldpk
    
    FI[1,1:n+3] <- - el*( 1- pkAk/eplkmo*lam*(1-eps) ) + tau * Tp *sum(Tl)
    FI[1:n+3,1] <- FI[1,1:n+3]
    
    ## d^2L/deps^2
    
    FI[2,2] <- el*sum(eplmoAl)/(1-eps)  + tau*sum(Te)^2
    
    ## d^2L/depsdpk
    
    FI[2,1:n+3] <- - el*lam/Ak  + tau*sum(Te) * Tp
    FI[1:n+3,2] <- FI[2,1:n+3]
    
    ## d^2L/dpk^2
    
    D <- (el*(1-eps)*lam^2/Ak/eplkmo + tau *  Tp^2)
    FI <- FI + diag(c(0,0,0,D ))
    
    ## use blockwise inverision
    FI <- N/(el-1)*FI
    
    D <- N/(el-1)*D
    D1 <- 1/D
    d <- -sum(D1)
    D2 <- c(-1,D1)
    D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
    
    A0 <-FI[1:2,1:2]
    B <- FI[1:2,3:(n+3)]
    A1 <- FI[1:n+3,1:n+3]
    DC <- D3 %*% t(B)
    A1 <- solve(A0 - (B %*% DC))
    
    iFI[1:2,1:2] <- A1
    C2 <- - DC %*% A1 
    iFI[0:n+3,1:2] <- C2
    B2 <- t(C2)
    iFI[1:2,0:n+3] <- B2
    iFI[0:n+3,0:n+3] <- D3 - DC %*% B2 
    iFI <- iFI[-3,-3]
  }
  iFI
}

lseq <- seq(0.1,2.5,0.1)
for(lam in lseq){
  adj <- exp(lam)*(exp(lam)-lam-1)/(exp(lam)-1)^2
  out <- c(out,adj^2*FI(N,lam,p,eps)[1,1])
  
}
FI(N,lam,p,eps)
out

x <-lseq/(1-exp(-lseq))

plot(x,sqrt(out)/x*100,type="l")


lseq <- seq(0.1,2.5,0.1)
out <-c()
for(lam in lseq){
  #adj <- exp(lam)*(exp(lam)-lam-1)/(exp(lam)-1)^2
  out <- c(out,FI(N,lam,p,eps)[3,3])
  
}

out

x <-lseq/(1-exp(-lseq))
print(x)
print(out)
plot(x,sqrt(out)/p[1]*100,type="l")


lseq <- seq(0.1,2.5,0.1)
out <-c()
for(lam in lseq){
  #adj <- exp(lam)*(exp(lam)-lam-1)/(exp(lam)-1)^2
  out <- c(out,FI(N,lam,p,eps)[2,2])
  
}

out

x <-lseq/(1-exp(-lseq))
print(x)
print(out)
plot(x,sqrt(out)/eps*100,type="l")
