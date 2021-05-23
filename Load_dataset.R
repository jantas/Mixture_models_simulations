
load_dataset <- function(dataset, gen.params){
  # dataset: iris, crabs, gen
  
  switch(dataset,
         
         "gen" = {
           
           n <- gen.params[[1]]
           K <- gen.params[[2]]
           p <- gen.params[[3]]
           AvgOmega <- gen.params[[4]] # Omega=pairwise overlap
           
           pdf.mixed.true <- MixSim(BarOmega = AvgOmega, p = p, K = K, sph = FALSE, hom = FALSE)
           # sph stands for "spherical"
           
           data <- simdataset(n, Pi = pdf.mixed.true$Pi, Mu = pdf.mixed.true$Mu, S = pdf.mixed.true$S)
           
           trueID <- data$id
           
           X <- data$X
           X <- matrix(X,byrow = F,ncol = p)
         },
         
         "iris" = {
           
           data <- iris
           n <- dim(data)[1]
           D <- dim(data)[2]-1
           
           tmp <- data$Species
           trueID <- rep(0,n)
           trueID[tmp=="setosa"] <- 1; trueID[tmp=="versicolor"] <- 2; trueID[tmp=="virginica"] <- 3;
           
           X <- data[,1:4]
           X <- matrix(as.numeric(unlist(X)),byrow = F,ncol = D)
         },
         
         "crabs" = {
           
           data <- crabs
           n <- dim(data)[1]
           
           trueID <- rep(0,n)
           trueID[data$sp=="B" & data$sex=="F"] <- 1;
           trueID[data$sp=="B" & data$sex=="M"] <- 2;
           trueID[data$sp=="O" & data$sex=="F"] <- 3;
           trueID[data$sp=="O" & data$sex=="M"] <- 4;
           
           X <- data[,4:8]
           D <- dim(X)[2]
           X <- matrix(as.numeric(unlist(X)),byrow = F,ncol = D)
         }
  )
  return(list(X, trueID))
}