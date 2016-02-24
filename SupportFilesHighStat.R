#Functions for multivariate smoothers
#Code based on: Smoothing with Mixed Model Software
#               BY LONG NGO AND M.P. WAND
###################################################
Mystd <- function(x) {(x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)}

default.knots <- function(x,num.knots)
{
   if (missing(num.knots))
      num.knots <- max(5,min(floor(length(unique(x))/4),35))
   return(quantile(unique(x),seq(0,1,length=
                  (num.knots+2))[-c(1,(num.knots+2))]))
}
####################################################
tps.cov <- function(r)
{     
    r <- as.matrix(r)
    num.row <- nrow(r)
    num.col <- ncol(r)
    r <- as.vector(r)
    nzi <- (1:length(r))[r!=0]
    ans <- rep(0,length(r))
    ans[nzi] <- r[nzi]^2*log(abs(r[nzi]))
    if (num.col>1) ans <- matrix(ans,num.row,num.col)
    return(ans)
}


####################################################
Matern.cov <- function(r)
{     
    r <- as.matrix(r)
    num.row <- nrow(r)
    num.col <- ncol(r)
    r <- as.vector(r)
    ans <- (1+abs(r)) * exp(-1 * abs(r) )    
    if (num.col>1) ans <- matrix(ans,num.row,num.col)
    return(ans)
}



GetXandZ <- function(x, y, K, KnotPositions){
  X <- cbind(rep(1, length(y)),x, y)
  #X <- cbind(x, y)
  #X <- rep(1, length(y))
 
  #Now the Z matrix 
  dist.mat <- matrix(0,K,K)
  dist.mat[lower.tri(dist.mat)] <- dist(KnotPositions)
  dist.mat <- dist.mat + t(dist.mat)  



  #Max dist rho
#  dist.matxy <- matrix(0,length(x),length(y))
#  dist.matxy[lower.tri(dist.matxy)] <- dist(cbind(x,y))
#  dist.matxy <- dist.matxy + t(dist.matxy)   
#  rho <- max(dist.matxy)

  Omega <- tps.cov(dist.mat)   #/rho)  #K by K
  
  #Omega <- Matern.cov(dist.mat/rho)  #K by K

  #Get the Omega ^ -1/2 term...
  svd.Omega <- svd(Omega)
  sqrt.Omega <- t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
  #Omega.minhalf <- solve(sqrt.Omega)

  #First term in Z
  diffs.1 <- outer(x, KnotPositions[,1],"-")
  diffs.2 <- outer(y, KnotPositions[,2],"-")
  dists <- sqrt(diffs.1^2 + diffs.2^2)   
  
  #Zr <- Matern.cov(dists/rho)  
 
  Zr <- tps.cov(dists)   #/rho)  
  Z <- t(solve(sqrt.Omega, t(Zr)))
    
  list(X=X, Z=Z)
}

