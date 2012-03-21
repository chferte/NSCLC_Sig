### brig demo normalization supervised workflow

source("http://sage.fhcrc.org/CRAN.R")
pkgInstall("snm")



library(metaGEO)
library(snm)
seed <- 1234
set.seed(seed)
np <- 25000
na <- 50
gmeans <- rchisq(np,1,2)
gmeans[gmeans>15] <- runif(sum(gmeans>15),15,16)
data <- matrix(gmeans,nr=np,nc=na)
bio.var <- data.frame(groups=rep(c("A","B"),each=25))
adj.var <- data.frame(batches=rep(c("A","B","C","D","E"),times=10),
                      age=rnorm(50,1,0.5))
int.var <- data.frame(array=factor(1:50))
group.effect <- snm:::sim.probe.specific(data, bio.var$groups, 0.1, list(func=rnorm,params=c(mean=1,sd=0.3)))
batches.effect <- snm:::sim.probe.specific(data, adj.var$batches, 0.3, list(func=rnorm,params=c(mean=0,sd=0.3)))
age.effect <- snm:::sim.probe.specific(data, adj.var$age, 0.2, list(func=rnorm, params=c(mean=1,sd=0.1)))
M <- data + group.effect + batches.effect + age.effect
array.effect <- snm:::sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=2)))
E <- matrix(rnorm(length(data),0,0.25), nr=nrow(data), nc=ncol(data))
Y <- M + array.effect + E
true.nulls <- which(group.effect[,1] == group.effect[,26])
plot2 <- function(x,y, ...) {
  plot(Y[,c(x,y)], ...);
  abline(0,1,col="red",lwd=3,lty=2)
}
png(file="samples_33_vs_44.png")
plot2(44,33,xlab="Sample 44", ylab="Sample 33")
dev.off()
u <- fs(Y)
X <- model.matrix(~bio.var$groups)
Z <- model.matrix(~adj.var$age + adj.var$batches)
int.var <- data.frame(array=factor(1:ncol(Y)))
snm.fit <- snm(Y, bio.var=X, adj.var=Z, int.var=int.var, rm.adj=TRUE)
ks.test(snm.fit$pval[true.nulls], "punif")
hist(snm.fit$pval[true.nulls])