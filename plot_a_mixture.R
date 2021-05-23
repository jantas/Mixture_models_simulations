rm(list=ls())
library(RColorBrewer)


n <- 5000
K <- 10
p <- 2

par(mfrow=c(1,3))

for (omg in c(0.01,0.05,0.10)) {
  set.seed(15)
  pdf.mixed.true <- MixSim(BarOmega = omg, p = 2, K = 10, sph = FALSE, hom = FALSE)
  # sph stands for "spherical"
  data <- simdataset(n, Pi = pdf.mixed.true$Pi, Mu = pdf.mixed.true$Mu, S = pdf.mixed.true$S)
  
  trueID <- data$id
  trueID <- factor(trueID)
  jColors <- data.frame(trueID = levels(trueID),color = I(brewer.pal(nlevels(trueID), name = 'Paired')))
  
  X <- data$X
  X <- matrix(X,byrow = F,ncol = p)
  
  plot(X[,1], X[,2], col=jColors$color[match(trueID, jColors$trueID)], xlim=c(-0.5,1.5), ylim=c(-0.25,1.5), xlab = expression('x'[1]), ylab =  expression('x'[2]))
}
