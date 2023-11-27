library(ggplot2)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#D55E00",
               "#A6761D", "#332288","#1f78b4","#66A61E", "#E6AB02","#999999","#33a02c")	

xax <- expression(frac(lambda,1-e^- lambda))

lt <- c("IDM","OM")

lambda <- seq(0.1, 2.5, 0.1)
eps <- c(0.015, 0.1, 0.15)
ssize <- sort(c(50,100,150,200), decreasing = TRUE)
linfreq_missing_data <- list(c(0.5, 0.5), c(0.75, 0.25),
                             c(1/3, 1/3, 1/3), c(0.75, 0.15, 0.1),
                             c(0.25, 0.25, 0.25, 0.25), c(0.75, 0.15, 0.05, 0.05),
                             c(0.2, 0.2, 0.2, 0.2, 0.2), c(0.75, 0.1, 0.05, 0.05, 0.05))

path<- "C:/Users/meraj.hashemi/Downloads/I/PhD/"


###################### reading data #########################

ReaddatanestedN <- function(n,p,path,folder){
  data <- read.table(paste(path,folder,"/data-n", n, "-maxfreq", p,".txt",sep="",collapse=""),sep="\t")
  dat <-lapply(strsplit(as.matrix(data),split=" "),as.double)
  dat
}

###################### top title #########################

toptitle <- function(pp,n){
  if (length(pp[duplicated(pp)]) > 0 && length(pp[!duplicated(pp)]) != 1){
    ppt<-rev(as.vector(table(sort(pp,decreasing=TRUE))))
    ppu<-unique(sort(pp,decreasing=TRUE))
    title<-""
    t<-1
    for (i in 1:(length(ppt))){
      a<-t:(t+ppt[i] - 1)
      tit <- paste(",\"p\"[",a,"],\"=","\"",sep="",collapse="")
      tit <- paste(substring(tit,3,nchar(tit)-1),ppu[i],"\", \"",", ",collapse="",sep="")
      title<-paste("paste(\"",substring(title,8,(nchar(title)-2)),substring(tit,1,nchar(tit)),"\"",")",collapse="",sep="")
      t<-t+length(a)
      if (i==length(ppt)){title<-paste("paste(\"",substring(title,8,(nchar(title)-5)),")",collapse="",sep="")}
    }
  }else if (length(pp[!duplicated(pp)]) == 1){
    a<-1:n
    title <- paste(",\" p\"[",a,"],\"=","\"",sep="",collapse="")
    title <- paste("paste(\"",substring(title,4,nchar(title)-1),pp[1],"\"",")",collapse="")
  }else{
    a<-1:n
    title <- paste(",\","," p\"[",a,"],\"=",pp,"\"",sep="",collapse="")
    title <- paste("paste(\"",substring(title,4,nchar(title)),")",collapse="",sep="")
  }
  title
}

################# plot configuration #########################

plot_beauitfy <- function(p, title, yay) {
  p <- p + guides(color = guide_legend(override.aes = list(size = 1))) 
  p <- p + guides(linetype = guide_legend(override.aes = list(colour = "black")))
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
  p <- p + guides(linetype = guide_legend(order = 1)) 
  p <- p + guides(color = guide_legend(override.aes = list(size = 1), order = 2)) 
  p <- p + guides(fill = guide_legend(override.aes = list(color = "transparent"), order = 3))
  p <- p + labs(title=parse(text=title))
  p <- p + labs(x=xax, y = yay)
  p <- p + scale_colour_manual(values=cbPalette) 
  p
}

#######################create plot output ##########################

create_out <- function(dd, lln, i, ld, n, ev, lfreq) {
  data <- dd[(lln*(i-1) + 1):(lln*i), ]
  psi_org <- data[,1]/(1-exp(-data[,1]))
  psi <- rep(psi_org,2)
  lam <- rep(data[,1],2)
  N <- rep(data[,ld - n],2)
  e <- data[,2]
  if (ev == 'Bias') {
    dat1 <- as.vector(data[,6]) 
    dat2 <- as.vector(data[,7])
    dat <- c(dat1, dat2)
    dat <- 100*(dat/psi - 1)
    yay <- "relative bias in % (average MOI)"
  } 
  else if (ev == 'cv') {
    dat1 <- data[,9 + 5*n + 4] 
    dat2 <- data[,9 + 5*n + 5] 
    dat <- c(dat1, dat2)
    dat <- 100*sqrt(dat)/psi
    yay <- "CV in % (average MOI)"
  }
  else if (ev == 'mse') {
    dat1_mean <- as.vector(data[,6]) 
    dat2_mean <- as.vector(data[,7])
    dat1_var <- data[,9 + 5*n + 4] 
    dat2_var <- data[,9 + 5*n + 5] 
    dat1 <- (dat1_mean - psi_org)^2 + dat1_var
    dat2 <- (dat2_mean - psi_org)^2 + dat2_var
    dat <- c(dat1, dat2)
    dat <- 100*sqrt(dat)/psi
    yay <- "NRMSE in % (average MOI)" 
  }
  else if (ev == 'euclidean') {
    p_true <- matrix(lfreq, length(psi_org), n, byrow = TRUE)
    dat1 <- data[,7 + (1:n)] 
    dat2 <- data[,7 + n + (1:n)]
    dateu1 <- sqrt(rowSums((dat1 - p_true)^2))
    dateu2 <- sqrt(rowSums((dat2 - p_true)^2))
    dat <- c(dateu1, dateu2)
    yay <- "Euclidean Norm (lineage frequencies)" 
  }
  else if (ev == 'freq-bias') {
    p_true <- rep(lfreq[1], length(psi_org))
    dat1 <- data[,7 + 1] 
    dat2 <- data[,7 + n + 1]
    dat <- c(dat1, dat2)
    dat <- 100*(dat/rep(p_true,2) - 1)
    yay <- "relative bias in % (dominant lineage)" 
  }
  else if (ev == 'freq-cv') {
    p_true <- rep(lfreq[1], length(psi_org))
    dat1 <- data[,9 + 5*n + 5 + 1] 
    dat2 <- data[,9 + 5*n + 5 + n + 1]
    dat <- c(dat1, dat2)
    dat <- 100*sqrt(dat)/rep(p_true,2)
    yay <- "CV in % (dominant lineage)" 
  }
  else if (ev == 'eps-bias') {
    dat <- data[,3]
    dat <- 100*(dat/e - 1)
    yay <- expression("relative bias in % ("~epsilon~")")
  }
  else if (ev == 'eps-cv') {
    dat <- data[,9 + 5*n + 1]
    dat <- 100*sqrt(dat)/e
    yay <- expression("CV in % ("~epsilon~")")
  }
  out <- cbind(N, psi, e, dat)
  list(out, yay)
}


#######################plot relative bias/cv #######################

plotMissingData<-function(folder,linfreq, ssize, lambda, eps, ev, path){
  lln <- length(lambda)*length(ssize)
  for (p in linfreq){
    n <- length(p)
    print(p)
    lfreq <- as.vector(p)
    pp <- round(p,digits=2)
    d <- ReaddatanestedN(n, 100*pp[1],path,folder)
    ld <- length(d[[1]])
    dd <- t(matrix(unlist(d),ld, lln*length(eps)))
    for (i in 1:length(eps)) {
      out_list <- create_out(dd, lln, i, ld, n, ev, lfreq)
      out <- out_list[[1]]
      yay <- out_list[[2]]
      title <- toptitle(pp,n)  
      if (substr(ev,1,3) != 'eps') {
        outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], em=out[,4], ll = rep(lt, each = nrow(out)/2))
      }
      else {
        outt <- data.frame(N=as.factor(out[,1]), psi=out[,2], em=out[,4])
      }
      
      pdfname <- paste(path,folder,"/plots/psi-",ev,"-n", n, "-mfreq",100*pp[1],"-eps",eps[i]*100,".pdf",sep="",collapse="")
      pdf(file=pdfname,height=5)
      print((outt))
      if (substr(ev,1,3) != 'eps') {
        p <- ggplot(data=outt,aes(x=psi,y=em,colour=N, linetype = ll, fill = as.factor(i)))
      }
      else {
        p <- ggplot(data=outt,aes(x=psi,y=em,colour=N, fill = as.factor(i)))
      }
      p <- p + geom_line()
      p <- p + geom_point(color = "transparent")
      p <- p + geom_hline(yintercept=0,  color = "transparent")
      p <- plot_beauitfy(p, title, yay)
      if (substr(ev,1,3) != 'eps') {
        p <- p + scale_linetype_manual(name = "", values = c("solid","dashed"))
      }
      p <- p + scale_fill_manual(name  = bquote(epsilon == .(100*(eps[i]))~"%"), values = 1, labels = NULL)
      p 
      print(p)
      dev.off()
    }
  }
}