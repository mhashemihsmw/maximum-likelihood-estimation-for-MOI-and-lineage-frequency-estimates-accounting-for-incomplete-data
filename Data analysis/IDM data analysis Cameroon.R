library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/CameroonMcCollum2008.xlsx",2)
data

for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
}

output <- array(,c(2*(ncol(data)-1),5))
ct <- 0
for(j in 2:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  ct <- ct + (rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)

  out <- MLE(n0+Np,Nk,n_0=n0)
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  output[2*j-3,1:5] <- c(n0+Np,n0,out$`average MOI`,lcp,ucp)
    
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  
  output[2*j-2,1:5] <- c(n0+Np,NA,out$`average MOI`,lcp,ucp)
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}
sum(ct == 0)
output <- data.frame(output)
colnames(output) <- c("N","n0","MOI","lo","up")
output$loci <- rep(toupper(colnames(data))[-1],each=2)
output$model <- rep(c("IDM","OM"),ncol(data)-1)
output
output$n0[is.na(output$n0)] <-"" 


cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
p <- ggplot(data=output) 
p <- p +  geom_point(aes(x=loci,y=MOI, group=model,color=model),size=4,position=position_dodge(width=0.4))
p <- p + geom_errorbar(aes(x=loci,ymin=lo,ymax=up, group=model,color=model),size=0.5,width=0.65,position=position_dodge(width=0.4))
p <- p + geom_text(aes(x=loci,y=up,label = n0), vjust = -0.4,hjust = 0.82,size=5) # add labels x out of y genes x (y) 
p <- p + theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key=element_blank(),
               panel.border = element_rect(colour='black',fill=NA),
               axis.text = element_text(size = rel(1.3),color='black'),
               axis.text.x = element_text(angle=0,vjust=0.6),
               axis.title = element_text(size = rel(1.3)),
               plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
               legend.text = element_text(size = rel(1.3)),
               legend.title = element_text(size = rel(1.3)),
               legend.position="bottom",
               legend.direction="horizontal")
p <- p + scale_color_manual(values=cbPalette,name="Model:") + ylim(1,1.8)
p <- p +  annotate("text", x = 6.5, y =  1.8, label = expression(italic(N) == 331), hjust = 0, vjust = 1,size=5)
p <- p + labs(x="Loci",y="Mean MOI")
p



filename <- "/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/cameroon_results"
pdf(paste(filename,".pdf",sep=""),height=4.5,width=6)
print(p)
dev.off()


####### Cameroon 


data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/CameroonMcCollum2008.xlsx",2)
data

for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
}
sel1 <- unlist(lapply(data$sample.name, function(x) is.element(strsplit(x,"/")[[1]][2],c("04","05"))))
data <- data[!sel1,]

output <- array(,c(2*(ncol(data)-1),5))
for(j in 2:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)
  
  out <- MLE(n0+Np,Nk,n_0=n0)
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  output[2*j-3,1:5] <- c(n0+Np,n0,out$`average MOI`,lcp,ucp)
  
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  
  output[2*j-2,1:5] <- c(n0+Np,NA,out$`average MOI`,lcp,ucp)
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}


output <- data.frame(output)
colnames(output) <- c("N","n0","MOI","lo","up")
output$loci <- rep(toupper(colnames(data))[-1],each=2)
str(output)
output$loci <- factor(output$loci, labels= c("C3","E1", "Fr9","Fr10", "O1", "Q4","J6", "K6",	"U5",	"U6",	"L4",	"L5",	"J3",	"L1"))


output$model <- rep(c("IDM","OM"),ncol(data)-1)
output



cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
p <- ggplot(data=output) 
p <- p +  geom_point(aes(x=loci,y=MOI, group=model,color=model),size=4,position=position_dodge(width=0.4))
p <- p + geom_errorbar(aes(x=loci,ymin=lo,ymax=up, group=model,color=model),size=0.5,width=0.65,position=position_dodge(width=0.4))
p <- p + geom_text(aes(x=loci,y=up,label = n0), vjust = -0.4,hjust = 0.82,size=5) # add labels x out of y genes x (y) 
p <- p + theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key=element_blank(),
               panel.border = element_rect(colour='black',fill=NA),
               axis.text = element_text(size = rel(1.3),color='black'),
               axis.text.x = element_text(angle=0,vjust=0.6),
               axis.title = element_text(size = rel(1.3)),
               plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
               legend.text = element_text(size = rel(1.3)),
               legend.title = element_text(size = rel(1.3)),
               legend.position="bottom",
               legend.direction="horizontal")
p <- p + scale_color_manual(values=cbPalette,name="Model:") + ylim(0.8,2.5)
p <- p +  annotate("text", x = 12.5, y =  2.5, label = expression(italic(N) == 166), hjust = 0, vjust = 1,size=5)
p <- p + labs(x="Loci",y="Mean MOI")
p

output[output$model=="OM",]


filename <- "/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/cameroon_results"
pdf(paste(filename,".pdf",sep=""),height=4.5,width=6)
print(p)
dev.off()


####### Kenya




data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/KenyaMcCollum.xlsx",2,detectDates = TRUE)
data


for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
  if(is.na(data[k,3])){
    data[k,3] <- data[k-1,3]
  }
}

data <- data[!is.na(data[,4]),]
data <- data[as.Date(data[,4]) - as.Date("1994-03-31")<0,]
data <- data[,-(2:4)]
data
length(unique(data$sample))
j=2
output <- array(,c(2*(ncol(data)-1),5))
for(j in 2:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)
  
  out <- MLE(n0+Np,Nk,n_0=n0)
  out1 <- out
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  output[2*j-3,1:5] <- c(n0+Np,n0,out$`average MOI`,lcp,ucp)
  
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  
  output[2*j-2,1:5] <- c(n0+Np,NA,out$`average MOI`,lcp,ucp)
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}

output <- data.frame(output)
colnames(output) <- c("N","n0","MOI","lo","up")
output$loci <- rep(toupper(colnames(data))[-(1)],each=2)
output$loci <- factor(output$loci, labels= c("J6", "K6",	"U5",	"U6","U7",	"L4",	"L5",	"J3",	"L1"))

str(output)
#output$loci <- factor(output$loci, labels= c("C3","E1", "Fr9","Fr10", "O1", "Q4","J6", "K6",	"U5",	"U6",	"L4",	"L5",	"J3",	"L1"))


output$model <- rep(c("IDM","OM"),ncol(data)-1)
output



cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
p <- ggplot(data=output) 
p <- p +  geom_point(aes(x=loci,y=MOI, group=model,color=model),size=4,position=position_dodge(width=0.4))
p <- p + geom_errorbar(aes(x=loci,ymin=lo,ymax=up, group=model,color=model),size=0.5,width=0.65,position=position_dodge(width=0.4))
p <- p + geom_text(aes(x=loci,y=up,label = n0), vjust = -0.4,hjust = 1.2,size=5) # add labels x out of y genes x (y) 
p <- p + theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key=element_blank(),
               panel.border = element_rect(colour='black',fill=NA),
               axis.text = element_text(size = rel(1.3),color='black'),
               axis.text.x = element_text(angle=0,vjust=0.6),
               axis.title = element_text(size = rel(1.3)),
               plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
               legend.text = element_text(size = rel(1.3)),
               legend.title = element_text(size = rel(1.3)),
               legend.position="bottom",
               legend.direction="horizontal")
p <- p + scale_color_manual(values=cbPalette,name="Model:") + ylim(0.8,2.5)
p <- p +  annotate("text", x = 8.5, y =  2.5, label = expression(italic(N) == 43), hjust = 0, vjust = 1,size=5)
p <- p + labs(x="Loci",y="Mean MOI")
p

output[output$model=="OM",]


filename <- "/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/kenya_results"
pdf(paste(filename,".pdf",sep=""),height=4.5,width=6)
print(p)
dev.off()

#######


data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/CameroonMcCollum2008.xlsx",2)
data

for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
}
sel1 <- unlist(lapply(data$sample.name, function(x) is.element(strsplit(x,"/")[[1]][2],c("04","05"))))
data <- data[!sel1,]

output <- array(,c(2*(ncol(data)-1),5))
for(j in 2:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)
  
  out <- MLE(n0+Np,Nk,n_0=n0)
  lcp <- out$`MOI parameter lambda` - qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  lcp <- lcp/(1-exp(-lcp))
  ucp <- out$`MOI parameter lambda` + qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  ucp <- ucp/(1-exp(-ucp))
  output[2*j-3,1:5] <- c(n0+Np,n0,out$`average MOI`,lcp,ucp)
  
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")
  
  lcp <- out$`MOI parameter lambda` - qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  lcp <- lcp/(1-exp(-lcp))
  ucp <- out$`MOI parameter lambda` + qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  ucp <- ucp/(1-exp(-ucp))
  
  output[2*j-2,1:5] <- c(n0+Np,NA,out$`average MOI`,lcp,ucp)
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}

output <- data.frame(output)
colnames(output) <- c("N","n0","MOI","lo","up")
output$loci <- rep(toupper(colnames(data))[-1],each=2)
output$model <- rep(c("IDM","OM"),ncol(data)-1)
output
output$n0[is.na(output$n0)] <-"" 


cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
p <- ggplot(data=output) 
p <- p +  geom_point(aes(x=loci,y=MOI, group=model,color=model),size=4,position=position_dodge(width=0.4))
p <- p + geom_errorbar(aes(x=loci,ymin=lo,ymax=up, group=model,color=model),size=0.5,width=0.65,position=position_dodge(width=0.4))
p <- p + geom_text(aes(x=loci,y=up,label = n0), vjust = -0.4,hjust = 0.82,size=5) # add labels x out of y genes x (y) 
p <- p + theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key=element_blank(),
               panel.border = element_rect(colour='black',fill=NA),
               axis.text = element_text(size = rel(1.3),color='black'),
               axis.text.x = element_text(angle=0,vjust=0.6),
               axis.title = element_text(size = rel(1.3)),
               plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
               legend.text = element_text(size = rel(1.3)),
               legend.title = element_text(size = rel(1.3)),
               legend.position="bottom",
               legend.direction="horizontal")
p <- p + scale_color_manual(values=cbPalette,name="Model:") + ylim(0.9,2.4)
p <- p +  annotate("text", x = 6.5, y =  2.4, label = expression(italic(N) == 166), hjust = 0, vjust = 1,size=5)
p <- p + labs(x="Loci",y="Mean MOI")
p

output


filename <- "/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/cameroon_results"
pdf(paste(filename,".pdf",sep=""),height=4.5,width=6)
print(p)
dev.off()


#######


data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/CameroonMcCollum2008.xlsx",2)
data

for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
}

sel1 <- unlist(lapply(data$sample.name, function(x) is.element(strsplit(x,"/")[[1]][2],c("04","05"))))
data <- data[!sel1,]
output <- array(,c(2*(ncol(data)-1),5))
for(j in 2:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)
  
  out <- MLE(n0+Np,Nk,n_0=n0)
  lcp <- out$`MOI parameter lambda` - qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  lcp <- lcp/(1-exp(-lcp))
  ucp <- out$`MOI parameter lambda` + qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  ucp <- ucp/(1-exp(-ucp))
  output[2*j-3,1:5] <- c(n0+Np,n0,out$`average MOI`,lcp,ucp)
  
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")

  lcp <- out$`MOI parameter lambda` - qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  lcp <- lcp/(1-exp(-lcp))
  ucp <- out$`MOI parameter lambda` + qnorm(0.975)*sqrt(out$`inverse Fisher information`[1,1])
  ucp <- ucp/(1-exp(-ucp))
  
  output[2*j-2,1:5] <- c(n0+Np,NA,out$`average MOI`,lcp,ucp)
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}

output <- data.frame(output)
colnames(output) <- c("N","n0","MOI","lo","up")
output$loci <- rep(toupper(colnames(data))[-1],each=2)
output$model <- rep(c("IDM","OM"),ncol(data)-1)
output


cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")

p <- ggplot(data=output) 
p <- p +  geom_point(aes(x=loci,y=MOI, group=model,color=model),size=4,position=position_dodge(width=0.4))

p <- p + geom_errorbar(aes(x=loci,ymin=lo,ymax=up, group=model,color=model),size=0.5,width=0.65,position=position_dodge(width=0.4))
p <- p + geom_text(aes(x=loci,y=up,label = n0), vjust = -0.4,hjust = 0.82,size=5) # add labels x out of y genes x (y) 
p <- p + theme(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key=element_blank(),
               panel.border = element_rect(colour='black',fill=NA),
               axis.text = element_text(size = rel(1.3),color='black'),
               axis.text.x = element_text(angle=0,vjust=0.6),
               axis.title = element_text(size = rel(1.3)),
               plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
               legend.text = element_text(size = rel(1.3)),
               legend.title = element_text(size = rel(1.3)),
               legend.position="bottom",
               legend.direction="horizontal")
p <- p + scale_color_manual(values=cbPalette,name="Model:") + ylim(1,2.4)
#p <- p + scale_x_discrete(labels=eval(parse(text=aa)))  # here is the trick to get the greek characters in the axes
p <- p + labs(x="Loci",y="Mean MOI")

p



filename <- "/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/cameroon_results"
pdf(paste(filename,".pdf",sep=""),height=4.5,width=6)
print(p)
dev.off()



data <- read.xlsx("/Users/kristanschneider/Library/CloudStorage/GoogleDrive-mathmalaria@gmail.com/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Meraj/Missing Data/Data analysis/KenyaMcCollum.xlsx",2,detectDates = TRUE)
data


for(k in 2:nrow(data)){
  if(is.na(data[k,1])){
    data[k,1] <- data[k-1,1]
  }
  if(is.na(data[k,3])){
    data[k,3] <- data[k-1,3]
  }
}

data <- data[!is.na(data[,4]),]
data <- data[as.Date(data[,4]) - as.Date("1994-03-31")<0,]

length(unique(data$sample))
for(j in 5:ncol(data)){
  alleles <- unique(na.omit(data[,j]))
  splist <- split(data[,j],data[,1])
  mat01 <- t(array(unlist(lapply(splist,function(x) is.element(alleles,x))),c(length(alleles),length(splist))))
  n0 <- sum(rowSums(mat01)==0)
  Np <- sum(rowSums(mat01)!=0)
  Nk <- colSums(mat01)
  
  out <- MLE(n0+Np,Nk,n_0=n0)
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  
  print(n0)
  print(c(out$`average MOI`,lcp,ucp))
  
  out <- MLE(n0+Np,Nk,n_0=n0,model="OM")
  lcp <- out$`average MOI` - qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  ucp <- out$`average MOI` + qnorm(0.975)*sqrt(out$`inverse Fisher information adjusted for average MOI`[1,1])
  
  print(c(out$`average MOI`,lcp,ucp))
  
  
}