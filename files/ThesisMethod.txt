# Umea University 
# Kandidatuppsatsen / Bachelor Thesis
# Statistik C2
# Gaussian Graphical Models

# Packages
library(MASS) 
library(huge)
library(matlib) 
library(Matrix)
library(matrixcalc)
library(bootnet)
library(gRim)
library(gRbase)
library(qgraph)
library(CVglasso)
# library(parallel)
# library(muvis) 


# Generation function
Sparse_Covariance_Matrix <- function(variables, non_zero_entries,seed){
  set.seed(seed)
  library(Matrix)
  library(matrixcalc)
  var <- variables
  non_zero <- non_zero_entries
  
  # Uniform Function
  Random_Function <- function(n){
    runif(n, min = -1, max = 1)
  }
  
  # rsparsematrix from package matrixcalc
  (S9 <- rsparsematrix(var, var, nnz = non_zero,
                       rand.x = Random_Function, 
                       symmetric=TRUE))
  S9 <- as.matrix(S9)
  #  diag(S9) =  1
  library(matlib)
  if(det(S9) != 0 && (is.positive.definite(S9, tol=1e-4)) == T){
    return(S9)
  }
  else{
    while((det(S9)) == 0 && (is.positive.definite(S9, tol=1e-4)) == F){
      #print(paste("non zero entries:", non_zero))
      non_zero <- non_zero + 1
      (S9 <- rsparsematrix(var, var, nnz = non_zero,
                           rand.x = Random_Function,
                           symmetric=TRUE))
      S9 <- as.matrix(S9)
      diag(S9) =  1
    }
    diag(S9) =  1
    multi_return <- function() {
      my_list <- list("CovarianceMatrix" = S9, "Non Zero Entries" = non_zero)
      return(my_list)
    }
    Results <- multi_return()
  }
  return(Results)
}

# Function for finding right seed to fulfill positive definite matrix 
Omega.function <- function(x,y,z){
  x <- x; y <- y; z <- z
  library(matrixcalc)
  while(
    (is.positive.definite(Sparse_Covariance_Matrix(x, y, z)$CovarianceMatrix,
                          tol=1e-6)) == F ){
    print(paste("Seed:", z))
    z <- z + 1
  }
  data <- Sparse_Covariance_Matrix(x, y, z)$CovarianceMatrix
  non_zero <- Sparse_Covariance_Matrix(x, y, z)$'Non Zero Entries'
  
  multi_return <- function() {
    my_list <- list("CovarianceMatrix" = data,'Non Zero Entries' = non_zero ,"Seed" = z)
    return(my_list)
  }
  Results <- multi_return()
}

# Generated Co-variance Matrices
mydata1 <- Omega.function(5,8,2) 
mydata1$CovarianceMatrix # 5x5 matrix with 6 non zero entries

mydata2 <- Omega.function(15,10,14) 
mydata2$CovarianceMatrix # 15x15 matrix with 10 non zero entries

mydata3 <- Omega.function(25,15,10)
mydata3$CovarianceMatrix # 25x25 matrix with 16? non zero entries


mydata4 <- Omega.function(100,50,15460)
mydata4$CovarianceMatrix




Omegas <- list("Omega 1" = mydata1$CovarianceMatrix, # Determinant is 0.4227915
               "Omega 2" = mydata2$CovarianceMatrix, # Determinant is 0.01259296
               "Omega 3" = mydata3$CovarianceMatrix, # Determinant is 0.002286404
               "Omega 4" = mydata4$CovarianceMatrix, # Determinant is 9.469159e-13

# Removing excess files to keep Environment clean
rm(mydata1,mydata2,mydata3,mydata4)
attach(Omegas)
#Rouge plots of Generated Co-variance Matrices
huge.plot(`Omega 1`)
huge.plot(`Omega 2`)
huge.plot(`Omega 3`)
huge.plot(`Omega 4`)



# Simulation of Graph using Covariance Matrices
covariance.simulation <- function(iterations, observations, Omega, seed){
  obs <- observations
  iter <- iterations
  Omega <- Omega
  nu <- Sys.time()
  set.seed(seed)
  ZERO.matrix <- diag(x = 1, nrow = length(Omega[,1]), ncol = length(Omega[,1]))
  mean.mu <- integer(length(Omega[,1]))
  
  Result.huge.ebic    <- ZERO.matrix
  Result.huge.stars   <- ZERO.matrix
  Result.bootnet.ebic <- ZERO.matrix
  Result.CVGlasso.CV  <- ZERO.matrix
  Result.CVGlasso.AIC <- ZERO.matrix
  Result.CVGlasso.BIC <- ZERO.matrix
  Result.bootnet.step <- ZERO.matrix
    
    
    
   Model.selection.ebic   <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.stars  <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.boot   <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.bic    <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.CV     <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.AIC    <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.BIC    <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
   Model.selection.step   <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))

  library(matlib)
  if(det(Omega) != 0){
    inverse <- inv(Omega)
    
    library(parallel)
    library(huge)
    library(MASS)
    
    i <- 1
    while (i <= iter) {
      #Model.selection.ric   <- matrix(0,nrow = length(Omega[,1]), ncol = length(Omega[,1]))
      print(paste("This is iteration =", i))
      # Generation of Data
      multivar.data <- mvrnorm(obs , mean.mu, Sigma = inverse, tol = 1e-6) # changed to e-14 from e-6
      multivar.df   <- as.data.frame(multivar.data)
      
      # # Model estimation with huge package
      GGM_Model     <- huge(multivar.data, method = "glasso", cov.output = F)
      # 
      # # Model Estimation with bootnet
      GGM_Model2    <- estimateNetwork(multivar.df,
                                       default = "EBICglasso",
                                       tuning = 0.5)
      GGM_Model3 <- estimateNetwork(multivar.df,default = c("ggmModSelect"),
                                    tuning = 0, start = c("full"), 
                                    considerPerStep = c("subset"),
                                    nCores = 4,
                                    stepwise = T)


      # # Optimal covariance selection with huge package, "ebic", "stars" and "ric"
      Model.selection.ebic  <- huge.select(GGM_Model, criterion = "ebic")$opt.icov
      Model.selection.stars <- huge.select(GGM_Model, criterion = "stars")$opt.icov
      Model.selection.ric   <- huge.select(GGM_Model, criterion = "ric")$opt.icov
       
      # # Optimal covariance selection with bootnet package using ebic
      Model.selection.boot  <- GGM_Model2$graph
      # 
      # # Optimal covariance selection with qgraph package using stepwise bic
      Model.selection.step <- GGM_Model3$graph
      # 
      Model.selection.CV   <- CVglasso(multivar.data, maxit = 10000,
                                       lam.min.ratio = 1e-2,
                                       nlam = 10,
                                       crit.cv = "loglik",
                                       cores = 1,
                                       K = 5)$Omega

      Model.selection.AIC  <- CVglasso(multivar.data, maxit = 10000,
                                       lam.min.ratio = 1e-2,
                                       nlam = 10,
                                       crit.cv = "AIC",
                                       cores = 1,
                                       K = 5)$Omega

      Model.selection.BIC  <- CVglasso(multivar.data, maxit = 10000,
                                       lam.min.ratio = 1e-2,
                                       nlam = 10,
                                       crit.cv = "BIC",
                                       cores = 1,
                                       K = 5)$Omega
    

      # Loops in order to categorize if edge is present or not. 
      # Element is put as 1 if edge is present, otherwise 0,
      # and is then summed over S iterations.
      
      # # for CVglasso
      for(w in 1:(length(Omega))){
        if(abs(Model.selection.CV[w]) > 0){ # threshold put to 0
          Model.selection.CV[w] = 1
        }else{
          Model.selection.CV[w] = 0
        }
      }
      
      
      for(d in 1:(length(Omega))){
        if(abs(Model.selection.AIC[d]) > 0){ # threshold put to 0
          Model.selection.AIC[d] = 1
        }else{
          Model.selection.AIC[d] = 0
        }
      }

      # 
      for(u in 1:(length(Omega))){
        if(abs(Model.selection.BIC[u]) > 0){ # threshold put to 0
          Model.selection.BIC[u] = 1
        }else{
          Model.selection.BIC[u] = 0
        }
      }
      
      # # For bootnet 
      for(k in 1:(length(Omega))){
        if(abs(Model.selection.boot[k]) > 0){ # threshold put to 0
          Model.selection.boot[k] = 1
        }else{
          Model.selection.boot[k] = 0
        }
      }
      
      # for bootnet qgraph stepwise 
      for(g in 1:(length(Omega))){
        if(abs(Model.selection.step[g]) > 0){ # threshold put to 0
          Model.selection.step[g] = 1
        }else{
          Model.selection.step[g] = 0
        }
      }
      
      # 
      # For huge, ebic, stars, ric
      for(p in 1:(length(Omega))){
        if(abs(Model.selection.ebic[p]) > 0){ # threshold put to 0
          Model.selection.ebic[p] = 1
        }else{
          Model.selection.ebic[p] = 0
        }
      }

      for(j in 1:(length(Omega))){
        if(abs(Model.selection.stars[j]) > 0){ # threshold put to 0
          Model.selection.stars[j] = 1
          }else{
            Model.selection.stars[j] = 0
          }
      }


      Result.huge.ebic       <- Model.selection.ebic  + Result.huge.ebic
      Result.huge.stars      <- Model.selection.stars + Result.huge.stars
      Result.bootnet.ebic    <- Model.selection.boot  + Result.bootnet.ebic
      Result.CVGlasso.CV     <- Model.selection.CV   + Result.CVGlasso.CV
      Result.CVGlasso.AIC    <- Model.selection.AIC  + Result.CVGlasso.AIC
      Result.CVGlasso.BIC    <- Model.selection.BIC  + Result.CVGlasso.BIC
      Result.bootnet.step    <- Model.selection.step + Result.bootnet.step

      
      multi_return <- function() {
        my_list <- list("Glasso CV Loglik"= Result.CVGlasso.CV,
                        "Glasso CV BIC"   = Result.CVGlasso.BIC,
                        "Glasso CV AIC"   = Result.CVGlasso.AIC,
                        "huge.EBIC"       = Result.huge.ebic,
                        "huge.STARS"      = Result.huge.stars,
                        "bootnet.EBIC"    = Result.bootnet.ebic,
                        "qgraph.stepwise" = Result.bootnet.step,
                        "Original Matrix" = Omega
                        )
        return(my_list)
      }
      Results <- multi_return()
      
      i <- i + 1
    }
    Time <- Sys.time() - nu
    
    multi_return2 <- function() {
      my_list <- list( "Results" = Results,
                       "Simulation Time"    = Time
      )
      return(my_list)
    }
    
    Results <- multi_return2()
    
    return(Results)
  }else{
    geterrmessage("Determinant is Zero, Select other input")
  }
}
Iterations   <- 1000
Observations <- 1000 # Or 5000
Simulated.Results1 <- covariance.simulation(Iterations, Observations, `Omega 1` , 1)
Simulated.Results1
Simulated.Results2 <- covariance.simulation(Iterations, Observations, `Omega 2` , 1)
Simulated.Results2 
Simulated.Results3 <- covariance.simulation(Iterations, Observations, `Omega 3` , 1)
Simulated.Results3 
Simulated.Results4 <- covariance.simulation(Iterations, Observations, `Omega 4` , 1)
Simulated.Results4



rm(Iterations, Observations)







##################################
##################################
#### Confusion Matrix Script Under
##################################
##################################




# Confusion Matrix Script
library(MASS)
library(class)
library(ISLR)
library(caret)

# Note ConfusionMatrix Algorithm only works for 1000 iterations, otherwise need revision.
ConfusionMatrix <- function(Simulation.Results, Original.Matrix, Simulation_Iterations){
  nu <- Sys.time()
  Iterations <- Simulation_Iterations
  # Converting all values to positive, to sort "Edge" or no "Edge"
  Original.Matrix <- abs(Original.Matrix)
  #  Simulation.Results <- abs(Simulation.Results)
  # Removing Edges in the middle
  diag(Original.Matrix)    <- 0
  diag(Simulation.Results) <- 0
  
  Contrast.Matrix <- rep(rep("Edge", (length(Original.Matrix))), Iterations)
  Contrast.Matrix[Original.Matrix == 0] <- "No Edge"
  
  # Contrasting the matrix
  Contrast.Matrix <- as.factor(Contrast.Matrix)
  contrasts(Contrast.Matrix)
  
  # Converting into Numeric
  Contrast.Matrix <- as.numeric(Contrast.Matrix)
  Contrast.Matrix[Contrast.Matrix == 2] <- 0
  dim(Contrast.Matrix) <- c(length(Original.Matrix), Iterations)
  Contrast.Matrix
  
  d <- NULL
  
  # Summing up all numerical entries
  for(i in 1:(length(Contrast.Matrix[,1]))){
    
    d[i] <- sum(Contrast.Matrix[i,])
  }
  
  dim(d) <- c(length(Original.Matrix[,1]), length(Original.Matrix[,1]))
  
  # Functions for TP,FP,FN,TN
  TruePositive  <- function(a,res){
    a <- a
    res <- res
    out <- NULL
    for(p in 1:(length(res))){
      if(a[p] == 0){
        out[p] <- 0
      }else{
        if(res[p] > a[p]){
          out[p] <- a[p]
        }else{
          out[p] <- res[p]
        }
      }
    }

    out <- (sum(out)/2) # Divide by 2, because of matrix symetry
    return(out)
  }
  FalsePositive <- function(d,res){
    d <- d
    res <- res
    out <- NULL
    for(p in 1:(length(res))){
      if(res[p] == 0){
        out[p] <- 0
      }else{
        if(d[p] == 0){
          out[p] <- res[p]
        }else{
          out[p] <- 0
        }
      }
    }

    out <- (sum(out)/2) # Divide by 2, because of matrix symetry
    return(out)
  }
  FalseNegative <- function(d,res){
    d <- d
    res <- res
    out <- NULL
    for(p in 1:(length(res))){
      if(d[p] == 0 && res[p] == 0){
        out[p] <- 0
      }else{
        if(d[p] != 0){
          out[p] <- ((d[p]) - (res[p]))
        }else{
          out[p] <- 0 # Might be an error
        }
      }
    }

    out <- (sum(out)/2) # Divide by 2, because of matrix symetry
    return(out)
  }
  TrueNegative  <- function(d,res,Iter){
    d <- d
    res <- res
    Iter <- Iter
    out <- NULL
    for(p in 1:(length(res))){
      if(d[p] == 0 && res[p] == 0){ # if both prediction and org. matrix = 0, put 1000 for iter
        out[p] <- Iter
      }else{
        if(d[p] == 0){
          out[p] <- (Iter - (res[p]))
        }else{
          out[p] <- 0
        }
      }
    }
    dim(out) <- c((length(res[,1])),(length(res[,1])))
    diag(out) <- 0
    out <- (sum(out)/2)
    return(out)
  }
  
  TP <- TruePositive(d,Simulation.Results)
  FP <- FalsePositive(d,Simulation.Results)
  FN <- FalseNegative(d,Simulation.Results)
  TN <- TrueNegative(d,Simulation.Results,Iterations)
  
  TPR <- (TP / (TP + FN))
  FPR <- (FP / (FP + TN))
  TNR <- (TN / (TN + FP))
  FNR <- (FN / (FN + TP))
  ACC <- ((TP + TN)/ ((TP + TN + FP + FN)))
  
  
  
  
  Time <- Sys.time() - nu
  
  multi_return <- function() {
    my_list <- list("Sensitivitet / TPR"  = TPR,
                    "Specificitet / TNR"  = TNR,
                    "Accuracy"            = ACC,
                    "False Negative Rate" = FNR,
                    "False Positive Rate" = FPR,
                    "TP" = TP,
                    "FP" = FP,
                    "TN" = TN,
                    "FN" = FN,
                    "Time" = Time
    )
    return(my_list)
  }
  Results <- multi_return()
  return(Results)
}



# Simulation 1
attach(Simulated.Results1$Results)
Conf.Results.Huge.EBIC      <- ConfusionMatrix(Simulation.Results    = huge.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.Huge.STARS    <- ConfusionMatrix(Simulation.Results     = huge.STARS,
                                              Original.Matrix        = `Original Matrix`,
                                              Simulation_Iterations  = 1000)
Conf.Results.bootnet.EBIC   <- ConfusionMatrix(Simulation.Results    = bootnet.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.logli <- ConfusionMatrix(Simulation.Results    = `Glasso CV Loglik`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.AIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV AIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.BIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV BIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.qgraph.step    <- ConfusionMatrix(Simulation.Results    = `qgraph.stepwise`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)


multi_return <- function() {
  my_list <- list("huge.EBIC"       = Conf.Results.Huge.EBIC,
                  "huge.STARS"      = Conf.Results.Huge.STARS,
                  "bootnet.EBIC"    = Conf.Results.bootnet.EBIC,
                  "CVglasso.Loglik" = Conf.Results.CVglasso.logli,
                  "CVglasso.AIC"    = Conf.Results.CVglasso.AIC,
                  "CVglasso.BIC"    = Conf.Results.CVglasso.BIC,
                  "qgraph.stepwise" = Conf.Results.qgraph.step
  )
  return(my_list)
}
Sim.Result1 <- multi_return()

# Remove excessive Results
rm(Conf.Results.Huge.EBIC,
   Conf.Results.Huge.STARS,
   Conf.Results.bootnet.EBIC,
   Conf.Results.CVglasso.logli,
   Conf.Results.CVglasso.AIC,
   Conf.Results.CVglasso.BIC,
   Conf.Results.qgraph.step 
   )

detach(Simulated.Results1$Results)

# Simulation 2
attach(Simulated.Results2$Results)
Conf.Results.Huge.EBIC      <- ConfusionMatrix(Simulation.Results    = huge.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.Huge.STARS     <- ConfusionMatrix(Simulation.Results    = huge.STARS,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.bootnet.EBIC   <- ConfusionMatrix(Simulation.Results    = bootnet.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.logli <- ConfusionMatrix(Simulation.Results    = `Glasso CV Loglik`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.AIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV AIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.BIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV BIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.qgraph.step    <- ConfusionMatrix(Simulation.Results    = `qgraph.stepwise`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
multi_return2 <- function() {
  my_list <- list("huge.EBIC"       = Conf.Results.Huge.EBIC,
                  "huge.STARS"      = Conf.Results.Huge.STARS,
                  "bootnet.EBIC"    = Conf.Results.bootnet.EBIC,
                  "CVglasso.Loglik" = Conf.Results.CVglasso.logli,
                  "CVglasso.AIC"    = Conf.Results.CVglasso.AIC,
                  "CVglasso.BIC"    = Conf.Results.CVglasso.BIC,
                  "qgraph.stepwise" = Conf.Results.qgraph.step
  )
  return(my_list)
}
Sim.Result2 <- multi_return2()

# Remove excessive Results
rm(Conf.Results.Huge.EBIC,
   Conf.Results.Huge.STARS,
   Conf.Results.bootnet.EBIC,
   Conf.Results.CVglasso.logli,
   Conf.Results.CVglasso.AIC,
   Conf.Results.CVglasso.BIC,
   Conf.Results.qgraph.step)

detach(Simulated.Results2$Results)

# Simulation 3
attach(Simulated.Results3$Results)
Conf.Results.Huge.EBIC      <- ConfusionMatrix(Simulation.Results    = huge.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.Huge.STARS     <- ConfusionMatrix(Simulation.Results    = huge.STARS,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.bootnet.EBIC   <- ConfusionMatrix(Simulation.Results    = bootnet.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.logli <- ConfusionMatrix(Simulation.Results    = `Glasso CV Loglik`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.AIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV AIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.BIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV BIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.qgraph.step    <- ConfusionMatrix(Simulation.Results    = qgraph.stepwise,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
multi_return3 <- function() {
  my_list <- list("huge.EBIC"       = Conf.Results.Huge.EBIC,
                  "huge.STARS"      = Conf.Results.Huge.STARS,
                  "bootnet.EBIC"    = Conf.Results.bootnet.EBIC,
                  "CVglasso.Loglik" = Conf.Results.CVglasso.logli,
                  "CVglasso.AIC"    = Conf.Results.CVglasso.AIC,
                  "CVglasso.BIC"    = Conf.Results.CVglasso.BIC,
                  "qgraph.stepwise" = Conf.Results.qgraph.step
  )
  return(my_list)
}
Sim.Result3 <- multi_return3()

# Remove excessive Results
rm(Conf.Results.Huge.EBIC,
   Conf.Results.Huge.STARS,
   Conf.Results.bootnet.EBIC,
   Conf.Results.CVglasso.logli,
   Conf.Results.CVglasso.AIC,
   Conf.Results.CVglasso.BIC,
   Conf.Results.qgraph.step)

detach(Simulated.Results3$Results)


# Simulation 4
attach(Simulated.Results4$Results)
Conf.Results.Huge.EBIC      <- ConfusionMatrix(Simulation.Results    = huge.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.Huge.STARS     <- ConfusionMatrix(Simulation.Results    = huge.STARS,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.bootnet.EBIC   <- ConfusionMatrix(Simulation.Results    = bootnet.EBIC,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.logli <- ConfusionMatrix(Simulation.Results    = `Glasso CV Loglik`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.AIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV AIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.CVglasso.BIC   <- ConfusionMatrix(Simulation.Results    = `Glasso CV BIC`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
Conf.Results.qgraph.step    <- ConfusionMatrix(Simulation.Results    = `qgraph.stepwise`,
                                               Original.Matrix       = `Original Matrix`,
                                               Simulation_Iterations = 1000)
multi_return4 <- function() {
  my_list <- list("huge.EBIC"       = Conf.Results.Huge.EBIC,
                  "huge.STARS"      = Conf.Results.Huge.STARS,
                  "bootnet.EBIC"    = Conf.Results.bootnet.EBIC,
                  "CVglasso.Loglik" = Conf.Results.CVglasso.logli,
                  "CVglasso.AIC"    = Conf.Results.CVglasso.AIC,
                  "CVglasso.BIC"    = Conf.Results.CVglasso.BIC,
                  "qgraph.stepwise" = Conf.Results.qgraph.step
  )
  return(my_list)
}
Sim.Result4 <- multi_return4()

# Remove excessive Results
rm(Conf.Results.Huge.EBIC,
   Conf.Results.Huge.STARS,
   Conf.Results.bootnet.EBIC,
   Conf.Results.CVglasso.logli,
   Conf.Results.CVglasso.AIC,
   Conf.Results.CVglasso.BIC,
   Conf.Results.qgraph.step)

detach(Simulated.Results4$Results)



multi_return5 <- function() {
  my_list <- list("Simulation 1"       = Sim.Result1,
                  "Simulation 2"       = Sim.Result2,
                  "Simulation 3"       = Sim.Result3,
                  "Simulation 4"       = Sim.Result4
  )
  return(my_list)
}
All.Results <- multi_return5()


rm(multi_return,
   multi_return2,
   multi_return3,
   multi_return4,
   multi_return5,
   Sim.Result1,
   Sim.Result2,
   Sim.Result3,
   Sim.Result4
)


All.Results










