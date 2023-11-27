p <- 1:2
p <- p/sum(p)
p <- c(0.5,0.45,0.05)
eps <- 0.001
lam <- 2
N <- 200
out <-c()
for(lam in seq(0.1,0.3,0.1)){
  n <-length(p)
  
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
  FI <- FI + diag(c(0,0,0,D ))# + tau * diag(c(0,0,0, Tp^2))
  
  FI <- N/(el-1)*FI
  FI[1:n+3,3] <- 1
  FI[3,1:n+3] <- FI[1:n+3,3] 
  # FI <- -FI
  
  #solve(FI[0:n+3,0:n+3])
  #det(FI)
  out <- c(out,sqrt(solve(FI)[1,1])*(1-eml)/lam)
  D <- N/(el-1)*D
  D1 <- 1/D
  d <- -sum(D1)
  D2 <- c(-1,D1)
  d
  D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
  

  #solve(FI)
  A0 <-FI[1:2,1:2]
  B <- FI[1:2,3:(n+3)]
  A1 <- FI[1:n+3,1:n+3]
  #FI
  #A0
  #B
  #FI[0:n+3,0:n+3]%*% solve(FI[0:n+3,0:n+3])
  solve(FI[0:n+3,0:n+3])
  
  DC <- D3 %*% t(B)
  #B1 <- B %*% D1
  A1 <- solve(A0 - (B %*% DC))
  
  iFI[1:2,1:2] <- A1
  C2 <- - DC %*% A1 
  iFI[0:n+3,1:2] <- C2
  B2 <- t(C2)
  iFI[1:2,0:n+3] <- B2
  iFI[0:n+3,0:n+3] <- D3 - DC %*% B2 
  iFI

  
  
}
iFI1 <- iFI
rowSums(iFI[1:n+3,1:n+3])
solve(FI)

#solve(FI)[-3,-3]
lseq <- seq(0.1,2,0.1)
plot(lseq/(1-exp(-lseq)),out,type="l")


#####
# eps =0
for(lam in seq(0.1,0.3,0.1)){
  n <-length(p)
  
  FI <- array(0,c(n+2,n+2))
  iFI <- array(0,c(n+2,n+2))
  el <- exp(lam)
  eml <-1/el
  elk <- exp(lam*p)
  
  eplkmo <- (elk-1)
  eplmoAl <- t(eplkmo/Ak)
  Tl <- eps * p * elk /Ak
  Te <- eplkmo/Ak
  Tp <- eps*lam*elk/Ak
  ## d^2L/dl^2
  FI[1,1] <- -1/(1-eml)+el*sum( p^2/eplkmo)
  

  ## d^2L/dldpk
  
  FI[1,1:n+2] <- - el*(1+ p/eplkmo*lam )
  FI[1:n+2,1] <- FI[1,1:n+2]
  
  ## d^2L/dpk^2
  
  D <- (el*lam^2/eplkmo)
  FI <- FI + diag(c(0,0,D ))
  
  FI <- N/(el-1)*FI
  FI[1:n+2,2] <- 1
  FI[2,1:n+2] <- FI[1:n+2,2] 
  
  # FI <- -FI
  
 
  out <- c(out,sqrt(solve(FI)[1,1])*(1-eml)/lam)
  D <- N/(el-1)*D
  D1 <- 1/D
  d <- -sum(D1)
  D2 <- c(-1,D1)
  d
  D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
  
  A0 <-FI[1,1]
  B <- FI[1,2:(n+2)]
  A1 <- FI[1:n+2,1:n+2]
  DC <- D3 %*% B
  A1 <- solve(A0 - (B %*% DC))
  
  iFI[1:2,1:2] <- A1
  C2 <- - DC %*% A1 
  iFI[0:n+2,1:1] <- C2
  B2 <- t(C2)
  iFI[1:2,0:n+2] <- B2
  iFI[0:n+2,0:n+2] <- D3 - DC %*% B2 
  iFI
  
  
  
}
iFI
iFI1

n <-length(p)

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


## use blockwise inverision
D <- (el*(1-eps)*lam^2/Ak/eplkmo + tau *  Tp^2)
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
iFI[-3,-3]

FI[3,1:n+3] <-1
FI[1:n+3,3] <-1
solve(FI+diag(c(0,0,0,D)))[-3,-3]



########## eps=0
FI <- array(0,c(n+2,n+2))
iFI <- array(0,c(n+2,n+2))
el <- exp(lam)
eml <-1/el
elk <- exp(lam*p)
eplkmo <- (elk-1)


## d^2L/dl^2
FI[1,1] <- -1/(1-eml)+el*sum( p^2/eplkmo)

## d^2L/dldpk

FI[1,1:n+2] <- - el*(1+ p/eplkmo*lam )
FI[1:n+2,1] <- FI[1,1:n+2]

## d^2L/dpk^2

D <- (el*lam^2/eplkmo)
FI <- FI + diag(c(0,0,D ))

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
A1 <- 1/(A0 - (B %*% DC))

iFI[1,1] <- A1
C2 <-  -DC %*% A1 
#C2 <-  -DC %*% A1 
iFI[0:n+2,1] <- C2
B2 <- t(C2)
iFI[1,0:n+2] <- B2
iFI[0:n+2,0:n+2] <- D3 - DC %*% B2 
iFI[-2,-2]
solve(FI)[-2,-2]



