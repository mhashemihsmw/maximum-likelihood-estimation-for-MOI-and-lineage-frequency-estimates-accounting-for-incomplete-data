
#-------------------------------------Missing Data--------------------------------------

########################################################################################
#-------------------------------------EM Algorithm--------------------------------------
########################################################################################

EM_MissingData <- function(N, N0, Nk, lambda_initial, eps_initial) {
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
    gg <- 1
    k <- 1
    ll0 <- 0
    ll <- 1
    while(abs(lamt - lamnext) +abs(epst - epsnext) + sqrt(sum((pkt - pnext)^2)) > thr1) {
        lamt <- lamnext
        epst <- epsnext
        pkt <- pnext
        nextiter <- EM_next_iteration(lamnext, pnext, epst, N, N0, Nk, Nnk, snk)
        pnext <- nextiter[[1]]
        epsnext<- nextiter[[2]] 
        wt <- nextiter[[3]]
        ntt <- wt/N
        lamnext <-  ntt + 1
        lamnext <- EM_lambda_Newton(lamnext, N, thr1, ntt)
        #ll0 <- ll
        #ll <- logl_MissingData(N, N0, Nk, lamnext, pnext, epsnext)
    }
    
    pnextz <- array(0,length(z))  
    pnextz[z > 0] <- pnext 
    psi <- lamnext/(1-exp(-lamnext))
    final <- list(epsnext, lamnext, psi, pnextz)
    final
}

########################################################################################

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

########################################################################################


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

########################################################################################


logl_MissingData <- function (N, N0, Nk, lame, pe, epse) {
    explp <- exp(lame*pe) - 1
    explpe <- epse*(explp) + 1
    l1 <- -N*log(exp(lame) - 1)
    l2 <- sum(Nk)*log(1-epse)
    l3 <- sum(Nk*log(explp))
    l4 <- sum((N - Nk)*log(explpe))
    l5 <- N0*log(1 - (prod(explpe))^(-1))
    ll <- l1 + l2 + l3 + l4 + l5
    ll
}


########################################################################################
#-------------------------------------------MLE-----------------------------------------
########################################################################################




MLECP_MissingData <- function(Nk,N,la) {
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

MLECP_MissingData(c(50,20,40),100)
########################################################################################
MLEP_MissingData <- function(Nk,N){
    
    lambdaE <- -sum(log(1-(Nk/N)))
    pE <- (-1/lambdaE)*log(1-Nk/N)
    
    c(lambdaE,pE)
}


########################################################################################
#-------------------------------------Missing Data Generator----------------------------
########################################################################################

RandomDataset_MissingData <- function(data, eps, N, n){
    ran <- runif(N*n)
    ran <- matrix((ran > eps)*1, N, n)
    ndata <- ran*data
    nkmiss <- colSums(ndata)
    nzero <- sum((rowSums(ndata) == 0)*1)
    list(nzero,nkmiss,data)
}


#--------------------------------------------------------------------------------------------------------------------
#
###for nested sample generator, generates sample of size N
innersamplegenerator_MissingData <- function(N0,Nk,N,data,e,lambda,p,n) {
    regular <- FALSE
    Nn <- N^n
    if (is.na(N0) == TRUE) {
        data <- mnom(N, lambda, p, n)
        result <- RandomDataset_MissingData(data, e, N, n)
        N0 <- result[[1]]
        Nk <- result[[2]]
        data <- result[[3]]
    }
    if (sum(Nk > 0) == 1) {
        regular <- FALSE
    }
    else if (N0 != 0 & (sum(Nk) == (N - N0) || abs(prod(1 - Nk/N) - N0/N) < 1/(Nn) ){  ### prod(1 - Nk/N) - N0/N == 0 can lead to numerical errors e.g., prod(1-c(30,30)/50)==8/50 gives FALSE, difference 2.775558e-17.   
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
        data <- mnom(N, lambda, p, n)
        result <- RandomDataset_MissingData(data, e, N, n)
        N0 <- result[[1]]
        Nk <- result[[2]]
        if (sum(Nk > 0) == 1) {
            regular <- FALSE
        }
        else if (N0 != 0 & (sum(Nk) == (N - N0) || abs(prod(1 - Nk/N) - N0/N) <1/(Nn)){
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

########################################################################################
#-------------------------------------Simulation--------------------------------------
########################################################################################

linfreq_missing_data <- list(c(0.5, 0.5), c(0.75, 0.25))

linfreq_missing_data <- list(c(1/3, 1/3, 1/3), c(0.75, 0.15, 0.1))

linfreq_missing_data <- list(c(0.25, 0.25, 0.25, 0.25), c(0.75, 0.15, 0.05, 0.05))

linfreq_missing_data <- list(c(0.2, 0.2, 0.2, 0.2, 0.2))

linfreq_missing_data <- list( c(0.75, 0.1, 0.05, 0.05, 0.05))


linfreq_missing_data <- list(c(0.5, 0.5), c(0.75, 0.25),
                             c(1/3, 1/3, 1/3), c(0.75, 0.15, 0.1),
                             c(0.25, 0.25, 0.25, 0.25), c(0.75, 0.15, 0.05, 0.05),
                             c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.75, 0.1, 0.05, 0.05, 0.05))



#simulation_MissingData(10, linfreq_missing_data, folder = '', path = )

simulation_MissingData <- function(S, linfreq, folder,path="C:/Users/meraj.hashemi/Downloads/"){
    #path <- "C:/Users/meraj.hashemi/Downloads/"
    s <- Sys.time()
    lambda <- seq(0.1, 2.5, 0.1)
    eps <- c(0.015, 0.1, 0.15)
    ssize <- sort(c(50,100,150,200), decreasing = TRUE)
    print(ssize)
    NN <- ssize[1]
    restN <- ssize[-1]
    sz <- length(ssize)
    for (p in linfreq) {
        print(p)
        pp <- 100*round(p[1],digits = 2)
        n <- length(p)
        for (e in eps) {
            print(e)
            for (lam in lambda) {
                print(lam)
                simfinal <- rep(list(array(NA, c(S, 7 + 5*n))), sz)
                sim <- 0
                #k <- 0
                # k<- rep(0,sz)
                dgen <- rep(list(list(0)), sz)
                while (sim < S) {
                    sim <- sim + 1
                    md <- innersamplegenerator_MissingData(NA, NA, NN, NA, e, lam, p, n)
                    N0 <- md[[1]]
                    Nk <- md[[2]]
                    data <- md[[3]]
                    em <- EM_MissingData(NN, N0, Nk, lam, e)
                    mle <-  EM_MissingData(NN-N0, 0, Nk, lam, e)
                    #mle <- MLECP_MissingData(Nk, NN - N0, lam)
                    
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
                        print(N)
                        t <- t + 1
                        dataN <- data[1:N,]
                        mdN <- RandomDataset_MissingData(dataN, e, N, n)
                        N0N <- mdN[[1]]
                        NkN <- mdN[[2]]
                        mdN <- innersamplegenerator_MissingData(N0N, NkN, N, dataN, e, lam, p, n)
                        N0N <- mdN[[1]]
                        NkN <- mdN[[2]]
                        print(c(N0N, NkN))
                        data <- md[[3]]
                        em <- EM_MissingData(N, N0N, NkN, lam, e)
                        mle <- MLECP_MissingData(NkN, N - N0N, lam)
                    
                        enem <- sqrt(sum((em[[4]]-p)^2))  ## Euclidian norm freqs - true parameter  - extended model
                        enom <- sqrt(sum((mle[[3]]-p)^2))  ## Euclidian norm freqs - true parameter - original model

                        elpk <- exp(-em[[2]] * em[[4]])
                        tpem <- (1-elpk)/(1-exp(-em[[2]]))                              ## true prevalence  - extended model
                        opem <- (1-em[[1]])*(1-elpk)/(1-prod(em[[1]]*(1-elpk)+elpk))    ## observed prevalence  -   extended model
                        tpom <- (1-exp(-mle[[1]] * mle[[3]]))/(1-exp(-mle[[1]]))        ## true prevalence  - original model

                        # ^eps, ^lam, ^lam*, ^psi, ^psi*, ^pp, ^pp* 
                        result <- c(em[[1]], em[[2]], mle[[1]], em[[3]], mle[[2]], em[[4]], mle[[3]],enem,enom,tpem,opem,tpom)
                        
                        #print(unlist(mle))
                        simfinal[[t]][sim,] <-  result
                    }
                }
               
                for(d in 1:sz){
                    final <- c(lam, e, apply(simfinal[[d]],2,mean), apply(simfinal[[d]],2,var), ssize[d], p)
                    write.table(t(final),paste(path,folder,"/data-n",toString(n), "-maxfreq",pp, ".txt",sep=""),
                                append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
                }
            }
        }
    }
    print(Sys.time() - s)
}





########################################################################################
#-----------------------------------likelihood plots------------------------------------
########################################################################################

plotf <- function (lam, N, ntt) {
    f <- lam - (1 - exp(-lam))*ntt
    print(f)
    plot(f ~ lam, type = "l")
    abline(h = 0, col = "red")
}

########################################################################################

plotdf <- function (lam, N, ntt,x) {
    f <- 1 - exp(-lam)*ntt
    print(f)
    plot(f ~ lam, type = "l")
    abline(v = x, col = "red")
}

########################################################################################

plot3L <- function (lam, eps, N, N0, Nk, p) {
    
    ll <- matrix(0,length(lam), length(eps))
    k <- 0
    for (l in lam) {
        k <- k + 1
        explp <- exp(l*p) - 1
        j <- 0
        for(e in eps) {
            j <- j + 1
            ll[k,j] <- -N*log(exp(l) - 1) + log(1 - e)*sum(Nk) + sum(Nk*log(explp)) + 
                sum((N - Nk)*log(e*explp + 1)) + N0*log(1 - prod(e*explp + 1)^(-1))
        }
    }
    pp <- plot_ly(y = lam, x = eps, z = ll) %>% add_surface(
        contours = list(
            z = list(
                show=TRUE,
                usecolormap=TRUE,
                highlightcolor="#ff0000",
                project=list(z=TRUE)
            ),
            start = -200,
            end= max(ll) + 0.5,
            size = 0.5
        )
    )
    
    a <- which(ll == max(ll),arr.ind = TRUE)
    print(max(ll), digits = 15)
    print(c(lam[a[1]], eps[a[2]]))
    pp
}

########################################################################################

plot2L <- function (lam, eps, N, N0, Nk, p) {
    
    ll <- matrix(0,length(lam), length(eps))
    k <- 0
    for (l in lam) {
        k <- k + 1
        explp <- exp(l*p) - 1
        j <- 0
        for(e in eps) {
            j <- j + 1
            ll[k,j] <- -N*log(exp(l) - 1) + log(1 - e)*sum(Nk) + sum(Nk*log(explp)) + 
                sum((N - Nk)*log(e*explp + 1)) + N0*log(1 - prod(e*explp + 1)^(-1))
        }
    }
    pp <- plotly::plot_ly(y = lam, x = eps, z = ll, type = "contour",
                          contours = list(start = min(ll) - 1,end= max(ll) + 0.5,size = 0.5))
    a <- which(ll == max(ll),arr.ind = TRUE)
    print(max(ll), digits = 15)
    print(c(lam[a[1]], eps[a[2]]))
    pp
}
########################################################################################

plotLlam <- function (lam, e, N, N0, Nk, p) {
    
    ll <- 1:length(lam)
    k <- 0
    for (l in lam) {
        k <- k + 1
        explp <- exp(l*p) - 1
        ll[k] <- -N*log(exp(l) - 1) + log(1 - e)*sum(Nk) + sum(Nk*log(explp)) + 
            sum((N - Nk)*log(e*explp + 1)) + N0*log(1 - prod(e*explp + 1)^(-1))
    }
    pp <- plotly::plot_ly(y = ll ,x = lam, mode = "lines", type = "scatter")
    print(ll)
    print(max(ll))
    print(which(ll == max(ll)))
    pp
}

########################################################################################

plotLeps <- function (l, eps, N, N0, Nk, p) {
    
    ll <- 1:length(eps)
    k <- 0
    for (e in eps) {
        k <- k + 1
        explp <- exp(l*p) - 1
        ll[k] <- -N*log(exp(l) - 1) + log(1 - e)*sum(Nk) + sum(Nk*log(explp)) + 
            sum((N - Nk)*log(e*explp + 1)) + N0*log(1 - prod(e*explp + 1)^(-1))
    }
    pp <- plotly::plot_ly(y = ll ,x = eps, mode = "lines", type = "scatter")
    print(ll)
    print(max(ll))
    a<- which(ll == max(ll))
    print(eps[a])
    pp
}













