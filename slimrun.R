
set.seed(33) #random number generator needs to be called ones
#run slim 
source("spatial_sim_w_pedigree_runner.R")


#create script
test <- spatial_sim_w_pedigree(ngen = 100, N = 20, output_file = "~/slimsim/out/test2.slim", )

#save slimscript to slim running folder
write(test, "~/slimrun/dummy.txt")

#call slim with the script
system("~/slimrun/slim ~/slimrun/dummy.txt")



#get output back into R
library(vcfR)
library(readr)


#read file and split into vcf per generation 
dummy <- read_lines("~/slimsim/out/test.slim")

start <- which(dummy=="#OUT")
end <- c(start[-1]-1,length(dummy))
for (i in 1:length(start))
{
  write_lines(dummy[(start[i]+8):end[i]],path=paste("~/slimsim/out/vcf",i,".vcf",sep=""))
}

files <- list.files("~/slimsim/out/", pattern=".vcf")
gls <- list()

for (i in 1:length(files)){
  dumvcf<- read.vcfR(file.path("~/slimsim/out",files[i]))
  gls[[i]]<- vcfR2genlight(dumvcf)
  names(gls)[i] <- paste("Gen_",dummy[start[i]+1], sep="")
  
  indnames <- as.numeric(unlist(strsplit(dummy[start[i]+4]," ")))
  xs <- as.numeric(unlist(strsplit(dummy[start[i]+2]," ")))
  ys <- as.numeric(unlist(strsplit(dummy[start[i]+3]," ")))
  parent12 <- matrix(as.numeric(unlist(strsplit(dummy[start[i]+6]," "))), ncol=2, byrow = T)
  parent12 <-  t(apply(parent12, MARGIN = 1, function(x) x[order(x)]))
  gp1234 <- matrix(as.numeric(unlist(strsplit(dummy[start[i]+7]," "))), ncol=4, byrow = T)
  gp1234 <- t(apply(gp1234, MARGIN = 1, function(x) x[order(x)]))
  
  gls[[i]]@other$ind.metrics <- data.frame(id=indnames, x=xs, y=ys, p1 =parent12[,1], p2 =parent12[,2], gp1 = gp1234[,1], gp2 = gp1234[,2], gp3 = gp1234[,3], gp4 = gp1234[,4])
  
  indNames(gls[[i]]) <- indnames
  pop(gls[[i]])<- rep("pop", nInd(gls[[i]]))
  gls[[i]]@other$xy <- data.frame(x=xs, y=ys)
  #put them on the equator
  gls[[i]]@other$latlong <- gls[[i]]@other$xy
  colnames(gls[[i]]@other$latlong) <- c("lat","lon")
}




par(mfrow=c(2,3))
lapply(gls, plot)



lapply(gls, function(x) plot(x@other$xy, asp=1, xlim=c(-1,1), ylim=c(-1,1)))



#distances between offsprings

#find offsprings

#full sibs
gg <- gls[[5]]

im <- gg@other$ind.metrics
sibs <- data.frame(sib1=NA, sib2=NA)
p12 <- paste(im$p1, im$p2)
cc<-1
 
  for (i in 1:length(p12))
  {
   
   nsib <- sum(p12 %in% p12[i]) -1
   if (nsib>0) {
     indi <- t(combn(which(p12 %in% p12[i]),2))
   for (ii in 1:nrow(indi))
     {
     sibs[cc,] = im$id[indi[ii,]]
     cc <- cc+1
     }
  }
  }

sibs <- t(apply(sibs,1, function(x) x[order(x)]))

sibs <- sibs[!duplicated(paste(sibs[,1], sibs[,2])),]
#all distances

sibs <- data.frame(sibs)
D <- as.matrix(dist(im[,c("x","y")]))
A <- gl.grm(gg)
diag(A)<- 0
sibs$dist=NA
sibs$glrm<- NA
#distances between 
for (i in 1:nrow(sibs))
{
  sibs$dist[i]<- D[im$id %in% sibs[i,1],im$id %in% sibs[i,2] ]
  sibs$glrm[i] <- A[im$id %in% sibs[i,1],im$id %in% sibs[i,2] ]
}

plot(density(as.dist(D)), ylim=c(0,4))
points(density(sibs$dist), type="l",col=4 )
summary(sibs$dist)
sigma <- sqrt(mean(sibs$dist^2)/2)
nrow(sibs)
sigma

plot(im[,c("x","y")], asp=1, xlim=c(-1,1), ylim=c(-1,1))
for (i in 1:nrow(sibs))
{
  xs<- c(im$x[im$id==sibs[i,1]], im$x[im$id==sibs[i,2]])
  ys<- c(im$y[im$id==sibs[i,1]], im$y[im$id==sibs[i,2]])
  lines(xs, ys, col="orange", lwd=2)
}

hist(A)
abline(v=mean(sibs$glrm))





