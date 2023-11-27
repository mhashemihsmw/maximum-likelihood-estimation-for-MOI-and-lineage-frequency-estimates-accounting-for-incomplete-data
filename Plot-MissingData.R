


'lambda <- seq(0.1, 2.1, 0.1)
eps <- c(0.01, 0.05, 0.1)
linfreq <- list(c(0.2,0.2,0.2,0.2,0.2))#, c(0.9,0.1), c(0.5,0.5), c(0.7,0.15,0.1,0.05))
ssize <- sort(c(40,50,70,100,150,200), decreasing = TRUE)'



############################################################

plotbiasMissData_psi<-function(folder,linfreq, ssize, lambda, eps){
  #path<-"C://Users/mhashemi/Google Drive/PHD/Papers/MissingData/simulations/"
  path<-"C://Users/Meraj/Google Drive/PHD/Papers/Missing Data/simulations/"
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <- length(p)
    p <- as.vector(p)
    pp <- round(p,digits=2)
    d <- ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    
    for (i in 1:length(eps)) {
      data <- dd[(lln*(i-1) + 1):(lln*i), ]
      psi <- data[,1]/(1-exp(-data[,1]))
      lam <- data[,1]
      N <- data[,ld - n]
      dat <- data[,c(6,7)] 
      dat <- 100*(dat/matrix(psi, length(psi),2) - 1)
      #dat <- 100*(dat/lam - 1)
      out <- cbind(N, psi, dat)
      #out <- cbind(N, lam, dat)
      xax <- expression(frac(lambda,1-e^- lambda))
      #xax <- expression(lambda)
      title <- toptitle(pp,n)
      outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], em=out[,3], mle = out[,4])
      pdfname <- paste(path,folder,"/plots/psi-Bias-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print(out)
      cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
                     "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	
      p <- ggplot(data=outt,aes(x=psi,y=em,group=N,colour=N))
      p <- p + geom_line()
      p <- p + geom_line(aes(y=mle),linetype= "dashed")
      p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
      p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
      p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
      p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
      p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
      p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
      p <- p + theme(axis.ticks = element_line(color = "black"))
      p <- p + theme(axis.title = element_text(size = rel(1.3)))
      p <- p + theme(legend.text = element_text(size = rel(1.3)))
      p <- p + theme(legend.title = element_text(size = rel(1.1)))
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      p <- p + labs(title=parse(text=title))
      p <- p + labs(x=xax,y="Bias in %")
      p <- p + scale_colour_manual(values=cbPalette)
      p 
      print(p)
      dev.off()
    }
  }
}



plotcvMissData_psi<-function(folder,linfreq, ssize, lambda, eps){
  #path<-"C://Users/mhashemi/Google Drive/PHD/Papers/MissingData/simulations/"
  path<-"C://Users/Meraj/Google Drive/PHD/Papers/Missing Data/simulations/"
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <- length(p)
    p <- as.vector(p)
    pp <- round(p,digits=2)
    d <- ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    
    for (i in 1:length(eps)) {
      data <- dd[(lln*(i-1) + 1):(lln*i), ]
      psi <- data[,1]/(1-exp(-data[,1]))
      lam <- data[,1]
      N <- data[,ld - n]
      dat <- data[,7 + 2*n + c(4,5)] 
      dat <- 100*sqrt(dat)/matrix(psi, length(psi),2) 
      #dat <- 100*(dat/lam - 1)
      out <- cbind(N, psi, dat)
      #out <- cbind(N, lam, dat)
      xax <- expression(frac(lambda,1-e^- lambda))
      #xax <- expression(lambda)
      title <- toptitle(pp,n)
      outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], em=out[,3], mle = out[,4])
      pdfname <- paste(path,folder,"/plots/psi-cv-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print(out)
      cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
                     "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	
      p <- ggplot(data=outt,aes(x=psi,y=em,group=N,colour=N))
      p <- p + geom_line()
      p <- p + geom_line(aes(y=mle),linetype= "dashed")
      p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
      p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
      p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
      p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
      p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
      p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
      p <- p + theme(axis.ticks = element_line(color = "black"))
      p <- p + theme(axis.title = element_text(size = rel(1.3)))
      p <- p + theme(legend.text = element_text(size = rel(1.3)))
      p <- p + theme(legend.title = element_text(size = rel(1.1)))
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      p <- p + labs(title=parse(text=title))
      p <- p + labs(x=xax,y="CV in %")
      p <- p + scale_colour_manual(values=cbPalette)
      p 
      print(p)
      dev.off()
    }
  }
}




plotbiasMissData_eps<-function(folder,linfreq, ssize, lambda, eps){
  #path<-"C://Users/mhashemi/Google Drive/PHD/Papers/MissingData/simulations/"
  path<-"C://Users/Meraj/Google Drive/PHD/Papers/Missing Data/simulations/"
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <- length(p)
    p <- as.vector(p)
    pp <- round(p,digits=2)
    d <- ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    
    for (i in 1:length(eps)) {
      data <- dd[(lln*(i-1) + 1):(lln*i), ]
      psi <- data[,1]/(1 - exp(-data[,1]))
      N <- data[,ld - n]
      dat <- data[,3] 
      dat <- 100*(dat/eps[i] - 1)
      out <- cbind(N, psi, dat)
      xax <- expression(frac(lambda,1-e^- lambda))
      title <- toptitle(pp,n)
      outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], emeps=out[,3])
      pdfname <- paste(path,folder,"/plots/eps-Bias-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print(out)
      cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
                     "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	
      p <- ggplot(data=outt,aes(x=psi,y=emeps,group=N,colour=N))
      p <- p + geom_line()
      p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
      p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
      p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
      p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
      p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
      p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
      p <- p + theme(axis.ticks = element_line(color = "black"))
      p <- p + theme(axis.title = element_text(size = rel(1.3)))
      p <- p + theme(legend.text = element_text(size = rel(1.3)))
      p <- p + theme(legend.title = element_text(size = rel(1.1)))
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      p <- p + labs(title=parse(text=title))
      p <- p + labs(x=xax,y="Bias in %")
      p <- p + scale_colour_manual(values=cbPalette)
      p 
      print(p)
      dev.off()
    }
  }
}




plotcvMissData_eps<-function(folder,linfreq, ssize, lambda, eps){
  #path<-"C://Users/mhashemi/Google Drive/PHD/Papers/MissingData/simulations/"
  path<-"C://Users/Meraj/Google Drive/PHD/Papers/Missing Data/simulations/"
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <- length(p)
    p <- as.vector(p)
    pp <- round(p,digits=2)
    d <- ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    
    for (i in 1:length(eps)) {
      data <- dd[(lln*(i-1) + 1):(lln*i), ]
      psi <- data[,1]/(1 - exp(-data[,1]))
      N <- data[,ld - n]
      dat <- data[,7 + 2*n + 1] 
      dat <- 100*sqrt(dat)/eps[i]
      out <- cbind(N, psi, dat)
      xax <- expression(frac(lambda,1-e^- lambda))
      title <- toptitle(pp,n)
      outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], emeps=out[,3])
      pdfname <- paste(path,folder,"/plots/eps-cv-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print(out)
      cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
                     "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	
      p <- ggplot(data=outt,aes(x=psi,y=emeps,group=N,colour=N))
      p <- p + geom_line()
      p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
      p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
      p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
      p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
      p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
      p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
      p <- p + theme(axis.ticks = element_line(color = "black"))
      p <- p + theme(axis.title = element_text(size = rel(1.3)))
      p <- p + theme(legend.text = element_text(size = rel(1.3)))
      p <- p + theme(legend.title = element_text(size = rel(1.1)))
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      p <- p + labs(title=parse(text=title))
      p <- p + labs(x=xax,y="CV")
      p <- p + scale_colour_manual(values=cbPalette)
      p 
      print(p)
      dev.off()
    }
  }
}





plotirregMissData<-function(folder,linfreq, ssize, lambda, eps){
  path<-"C://Users/Meraj/Google Drive/PHD/Papers/Missing Data/simulations/"
  #path<-"C://Users/Meraj/Google Drive/PHD/Papers/"
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <-length(p)
    p<-as.vector(p)
    pp<-round(p,digits=2)
    d<-ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    
    for (i in 1:length(eps)) {
      data <- dd[(lln*(i-1) + 1):(lln*i), ]
      lam <- data[,1]
      N <- data[,ld - n]
      dat <- data[,ld - n - 1] 
      dat <- dat/100
      out <- cbind(N, lam, dat)
      xax <- expression(frac(lambda,1-e^- lambda))
      title <- toptitle(pp,n)
      outt <- data.frame(N=as.factor(out[,1]), lam=out[,2], ir=out[,3])
      pdfname <- paste(path,folder,"/plots/irreg-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print(outt)
      cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
                     "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	
      p <- ggplot(data=outt,aes(x=lam,y=ir,group=N,colour=N))
      p <- p + geom_line()
      p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
      p <- p + theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
      p <- p + theme(panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank())
      p <- p + theme(panel.background = element_rect(colour='black',fill='white'))
      p <- p + theme(legend.key=element_rect(colour='white', fill='white'))
      p <- p + theme(axis.text = element_text(colour='black',size = rel(1.4)),title = element_text(size = rel(1.3)))
      p <- p + theme(axis.ticks = element_line(color = "black"))
      p <- p + theme(axis.title = element_text(size = rel(1.3)))
      p <- p + theme(legend.text = element_text(size = rel(1.3)))
      p <- p + theme(legend.title = element_text(size = rel(1.1)))
      p <- p + theme(plot.title = element_text(hjust = 0.5))
      p <- p + labs(title=parse(text=title))
      p <- p + labs(x=expression(lambda),y="ratio of irregular data in %")
      p <- p + scale_colour_manual(values=cbPalette)
      p 
      print(p)
      dev.off()
    }
  }
}