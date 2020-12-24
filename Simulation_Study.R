# Weerhandi et al. simulations
# Dec 11, 2020

# install.packages("parallel")
# install.packages("nlme")
library(parallel) # for mclapply()
library(nlme)

# Simulation Study
# Fix VarE at 10 and vary VarA to take values such as 1, 5, 10 
# Note: If factpr varoance os large no use of BLUP
# Keep sample size n from each gropu eaul (nothing much we learn by varying them)
# ng is number of groups; in trypical application it tend to be around 5 but could be as small as 3
simulate <- function(n = 10, ng = 3, VarA = .1, VarE = 1){
  N     <- n * ng # Total sample size
  beta  <- 10 
  alpha <- 1 # coeffieceint of X is not really important
  mu    <- 1 # common mean of Z is also not imporatnt to cxahnge
  X     <- rnorm(N, beta, 1) # Do not change this during simulation
  Z     <- rnorm(N, mu, .5)
  # Test Below showed that BLUP exacly agree with the one we used in one-way ANOVA case
  #if (M == 1) { 
  #  X <- rep(1, N)
  #  print("Testing the ANOVA case")
  #  Z <- rep(1, N)
  #}
  Segs <- rep(c(1 : ng), each = n)	
  # Siulate e and u and compare BLUP with eBLUPs
  outAll <- data.frame(matrix(0, ng, 4))	
  
  e <- rnorm(N, 0, sqrt(VarE))
  u <- rnorm(ng, mu, sqrt(VarA)) # u vaues around group mean
  uvec <- rep(u, each = n)
  Y <- alpha + beta * X + Z * uvec + e 
  lmx <- lm(Y ~ X, singular.ok = T) # Regular regression as in ANOVA
  del <- lmx$residuals # This is Y-betahat X
  # Also confirmed that it does not agrees with uhat=Yi_bar
  zprimedelta <- Z * del
  zz <- Z * Z # This is before aggregating by groups below; so not using Z'Z
  zprimedel.seg <- aggregate(zprimedelta, by = list(Segs), FUN = sum, simplify = TRUE)[,2]
  zprimez.seg <- aggregate(zz, by = list(Segs), FUN = sum, simplify = TRUE)[,2]
  #if (M == 1){ # Testing
  #    grpYbar <- aggregate(Y, by = list(Segs), FUN = mean) # Group means in ANOVA case before centering
  #    grpEst <- zprimedel.seg / zprimez.seg
  #    return(cbind(grpYbar[,2] - mean(Y), grpEst)) 
  #  }
  BLUP <- (zprimedel.seg / zprimez.seg) * (1 - VarE / (VarE + VarA * zprimez.seg))# where is the mu_hat
  dat <- cbind.data.frame(Y, X, Z, Segs)
  otmp <- getBLUP(Y, X, Z, as.factor(Segs), dat) 
  #r1 <- ng * (m - 1) + 1
  #r2 <- ng * m
  outAll[1 : ng, ] <- cbind.data.frame(BLUP, otmp)
  # Return all output and print MSEs
  #diff <- outAll[, 2 : 4] - outAll[, 1]
  #MSE <- colMeans(diff ^ 2)
  #names(MSE) <- c("REML", "ML", "Proposed")
  #print(MSE)
  names(outAll) <- c("BLUP", "REML", "ML", "Proposed")
  return(outAll)
  #return(MSE)
}

# run as out=getBLUP(Y,X,Z,Segement,data) or Default as getBLUP()
getBLUP <- function(Y, X, Z, segment, data){
  N <- length(Y)
  segs <- unique(segment)
  if (is.matrix(X)) {
    kx <- ncol(X)
  } else {
    kx <- 2
  }
  kx1 <- kx + 1 # number of X columns+2
  ng <- length(segs)
  #************************* First estimate rho values ******************
  lmxz <- lm(Y ~ X + Z : segment, singular.ok = T) # regular regression when 
  ilmxz <- as.matrix(summary(lmxz)$cov.unscaled)
  ni <- ncol(ilmxz)
  XX.inv <- ilmxz[kx1 : ni, kx1 : ni] # TO DO: Need to generalize for any Z
  # Now run GLSE to estimate beta and then u
  Zmat <- matrix(0, N, ng)
  for (j in 1 : ng) Zmat[is.element(segment, segs[j]), j] <- 1
  Zmat <- Zmat * Z
  X1 <- cbind(1, X)
  betx <- as.vector(lmxz$coef[1 : kx])
  # Use reduction of error going from lmx model to lmxz as in Sam's paper
  lmx <- lm(Y ~ X) # Regular regression when 
  err <- lmx$residuals
  # Just test below to show that residual from lmxz yield non-shrinked uhat0
  zprimez <- Z * Z ## This is z'z for all groups piled in same vector
  zprimedelta <- as.vector(Z * err) # Or do this using lmxz regression and (y-Xbeta_hat)
  # Both methods above yield similar etimates 
  zprimeresid.seg <- aggregate(zprimedelta, by = list(segment), FUN = sum, simplify = TRUE)[,2]
  zprimez.seg <- aggregate(zprimez, by = list(segment), FUN = sum, simplify = TRUE)[,2]
  uhat0 <- zprimeresid.seg / zprimez.seg 
  #****** Use uhat from lmxz regression, beacuse rho estimation part does not involve BLUP
  # Old test: uhat0=lmxz$coef[-c(1:kx)] # Note: err=Y-X1%*%betx applied to above formula also yiled same uhat0	
  muhat <- mean(uhat0)
  uhat <- uhat0 - muhat
  a <- length(table(segment)) - 1
  e <- lmxz$df.resid # Error df
  SSE <- sum(lmxz$residual ^ 2)
  mse <- SSE / e
  SSA <- anova(lmxz)[2, 2] # Noted that this is the same as SSA from lmx regression
  msa <- SSA / a
  M <- 1000 # Sample size for MonteCarlo Integreation
  Wae <- rf(M, a, e)
  Wae <- (e / (e - 2)) * Wae / mean(Wae)  
  good <- Wae >= (mse / msa)
  Wae <- Wae[good]
  cd <- a * (SSE / e) * mean(Wae)
  
  # Method without Eigen value decomposition: Just matrix inverse; This is less confusing and can do for any number of variance components
  rhohat <- getrho(cd, uhat, XX.inv)
  
  #**** Now compute eBLUP usineg Henderson Formulas*****************
  betahat <- getBetahat(Y, X1, Zmat, rhohat)
  e <- (Y - X1 %*% betahat)
  eBLUP  <- getUhat(Zmat, e, rhohat) # Return eBLUPs by REML and Proposed 
  eBLUP0 <- lme(Y ~ X, data = data, random = ~ 0 + Z | Segs)$coef$rand$Seg # REML
  eBLUP1 <- lme(Y ~ X, data, random = ~ 0 + Z | Segs, method = "ML")$coef$rand$Seg # ML
  eBLUPs <- cbind.data.frame(eBLUP0, eBLUP1, eBLUP)
  names(eBLUPs) <- c("REML", "ML", "Proposed")
  return(eBLUPs)
}

# Functions to get betahat and uhat using Henderson formulas
# TO DO: Allow tho to be matrix with few variance component
getBetahat <- function(Y, X, Z, rho){
  if (is.matrix(X)) N <- nrow(X)
  else N <- length(X)
  G <- diag(1, N) + rho * (Z %*% t(Z)) # Covarinace matrix of Y
  X1 <- X
  GI <- solve(G)
  betahat <- solve((t(X1) %*% GI %*% X1)) %*% (t(X1) %*% GI %*% Y) # This is correct and exactly te same as lme reuslts
  return(as.vector(betahat))
}
getUhat <- function(Z, e, rho){
  N <- nrow(Z)
  ng <- ncol(Z)
  rhoI <- 1 / rho # This need to be genarilzed as inverse
  uhat <- solve(t(Z) %*% Z + rhoI * diag(1, ng)) %*% t(Z) %*% e
  return(uhat)
}

# minimize function to estimate rho 
getrho <- function(cd, uhat, XX.inv){
  ng <- length(uhat)
  mindiff <- function(rho){
    XXInv <- solve(XX.inv + rho * diag(1, ng))
    return((t(uhat) %*% XXInv %*% uhat - cd) ^ 2)		
  }
  gout <- optimize(mindiff, lower = 0, upper = 50)
  return(gout$min)
}


# main simulation function to calculate MSEs
sim_par <- function(n = 10, ng = 3, VarA = .1, VarE = 1, rept = 10){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(11122023)
  res_list <- mclapply(1:rept, function(x){
    simulate(n, ng, VarA, VarE)
    }, mc.cores = detectCores() - 2)
  res <- do.call(rbind, res_list)
  #colnames(res) <- c("REML", "ML", "GM")
  colnames(res) <- c("BLUP", "REML", "ML", "GM")
  res
  }

# create table of the results
table_sim <- function(n = c(10, 50 ,100), 
                      ng = c(3, 5, 10), 
                      VarA = c(.1, 1, 10), 
                      VarE = c(1, 5, 10),
                      rept = 10000){
  
  #graph_mat <- matrix(0, length(n) * length(ng) * length(VarA) * length(VarE), 7)
  graph_mat <- matrix(0, length(n) * length(ng) * length(VarA) * length(VarE), 8)
  m <- 1
  for(i in 1:length(n)){
    for(j in 1:length(ng)){
      for(k in 1:length(VarA)){
        for(l in 1:length(VarE)){
          graph_mat[m,] <- c(n[i], ng[j], VarA[k], VarE[l], 
                             colMeans(sim_par(n[i], ng[j], VarA[k], VarE[l], rept)))
          m <- m + 1
        }
      }
    }
  }
  graph_mat <- data.frame(graph_mat)
  #colnames(graph_mat) <- c("n", "ng", "VarA", "VarE", "REML", "ML", "GM")
  colnames(graph_mat) <- c("n", "ng", "VarA", "VarE", "BLUP", "REML", "ML", "GM")
  #graph_mat_melted <- cbind(do.call("rbind", replicate(3, select(graph_mat, n, ng, VarA), simplify = FALSE)),
  #                          melt(select(graph_mat, REML, ML, GM)))
  # graph_mat_melted %>%
  #   ggplot(aes(x = factor(VarA), 
  #              y = value, 
  #              fill = variable)) +
  #   geom_boxplot(outlier.shape = NA) +
  #   ylim(c(0, 6)) +
  #   facet_grid(vars(n), vars(ng))
  
  #res_mat <- graph_mat_melted %>%
  #             group_by(n, ng, VarA, variable) %>%
  #             summarise(mean(value))
  #print(res_mat, n = dim(res_mat)[1])
  graph_mat
}
start <- Sys.time()
table_sim <- table_sim()
write.csv(table_sim, "table_sim.csv")
Sys.time() - start

