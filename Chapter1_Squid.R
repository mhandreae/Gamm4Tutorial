# R code for Chapter 1 in: Beginner's Guide to GAMM with R (2014)
# Zuur AF, Saveliev, AA, Ieno EN
# www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
##################################################################################

# Note:
# The code in this solution file was modified in May 2015. We originally (2012-2014) 
# used multivariate normal priors for regression parameters and random effects. 
# This is recommended by some. However, we noticed that mixing is not always
# good with multivariate normal priors, and we therefore changed the R code
# to univariate normal priors. The text in the chapter was modified accordingly.
# If you bought an e-copy of the GAMM book in June 2015 or later, then you will 
# see the updated text. If you bought a copy of the book in 2013, 2014, or 2015
# then you will see the old text (but there are only minor changes). The second
# print will contain updated text (and will be available from 2016 onwards).

 # Alain Zuur
 # May 2015, Newburgh, UK
 ########################
 





#Set the working directory
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/Chapter1")

#Load the data
Squid <- read.table(file = "SquidNorway.txt",
                    header = TRUE,
                    dec = ".")

#Check what we have
names(Squid)
str(Squid)
##################################################################################



##################################################################################
#Load packages and source files
library(lattice)
library(mgcv)
library(R2jags)

source("/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV6.R")
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/SupportFilesHighStat.R")
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstat.R")
##################################################################################



##################################################################################
#Figure 1.1. The easy way to make this graph
#You first need to source the Highstatlib.R file.
MyVar <- c("Lat", "Depth", "ML")
Myxyplot(Squid,MyVar,"d15N")





#And the more difficult way to make this graph:
Y    <- Squid$d15N
MyX  <- c("Lat", "Depth", "ML")
X    <- Squid[, MyX]

AllX  <- as.vector(as.matrix(Squid[,MyX]))
AllY  <- rep(Squid[,"d15N"] , length(MyX))
AllID <- rep(MyX, each = nrow(Squid))
    
xyplot(AllY ~ AllX|factor(AllID), col = 1,layout=c(3,1),
              xlab = list(label="Explanatory variables",cex = 1.5),
              ylab = list(label = "d15N", cex = 1.5),
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = T,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = 1)
                panel.loess(x, y, col = 1, lwd = 2)})
                
##################################################################################

 
  
##################################################################################
#Figure 1.2

M1 <- lm(d15N ~ Lat + ML, data = Squid)
E1 <- rstandard(M1)
par(mar = c(5,5,2,2))

plot(x = Squid$ML,
     y = E1,
     xlab = "ML",
     ylab = "Standardized residuals",
     cex.lab = 1.5, pch = 16)
abline(h = 0)

L1 <- loess(E1~Squid$ML)
I1 <- order(Squid$ML)
lines(x = sort(Squid$ML), 
      y = fitted(L1)[I1], 
      lwd = 3)

M2A <- gam(d15N ~ Lat + s(ML), data = Squid)
M2B <- gam(d15N ~ Lat + s(ML, bs = "cr"), data = Squid)
M2C <- gam(d15N ~ Lat + s(ML, bs = "ps"), data = Squid)
M2D <- gam(d15N ~ Lat + s(ML, bs = "cs"), data = Squid)
#################################################################




#################################################################
#Setion 1.2

Squid$Lat.std <- Mystd(Squid$Lat)
Squid$ML.std  <- Mystd(Squid$ML)

M3 <- gam(d15N ~ Lat.std + s(ML.std, bs = "cr"), data = Squid)
#Figure 1.3


par(mar = c(5,5,2,2))
plot(M3, resid = TRUE, pch = 16, cex = 0.7, cex.lab = 1.5)
summary(M3)


##################################################
#Section 1.3
#Linear spline
rhs <- function(x, TH) ifelse(x >= TH, x - TH, 0)
M4  <- lm(d15N ~ ML.std  + rhs(ML.std, 0.8) + Lat.std,
          data = Squid)
summary(M4)


#Figure 1.4
range(Squid$ML.std)
NewData <- data.frame(ML.std = seq(from = -1.48, to = 3.03, length = 100),
                      Lat.std = 0)
P4  <- predict(M4, newdata = NewData)


par(mar = c(5,5,3,3))
plot(x = Squid$ML.std,
     y = Squid$d15N, 
     pch = 16, 
     xlab = "Standardized ML",
     ylab = "d15N",
     cex.lab = 1.5)
abline(v = 0.83 , 
       lty = 2, 
       lwd = 2)     
lines(x = NewData$ML.std, 
      y = P4, 
      lwd = 5)


#Knots
probs <- seq(0, 1, length = 5)
QD    <- quantile(Squid$ML.std, probs)
QD

M5  <- lm(d15N ~ Lat.std  + ML.std  + 
                 rhs(ML.std, -0.7121410) +  
                 rhs(ML.std, -0.1667513) +
                 rhs(ML.std, 0.6419299) ,
          data = Squid)

##################################################################
#Figure 1.5
NewData <- data.frame(ML.std = seq(from = -1.48, to = 3.03, length = 100),
                      Lat.std = 0)
P5  <- predict(M5, newdata = NewData)



par(mar = c(5,5,3,3))
plot(x = Squid$ML.std,
     y = Squid$d15N, 
     pch = 16, 
     xlab = "Standardized ML",
     ylab = "d15N",
     cex.lab = 1.5)
lines(x = NewData$ML.std, 
      y = P5, 
      lwd = 5)
abline(v = -0.7121410, lty = 2, lwd = 2)     
abline(v = -0.1667513, lty = 2, lwd = 2)     
abline(v = 0.6419299, lty = 2, lwd = 2) 

#############################################################




#############################################################     
#Section 1.4
#Now fit the same model in JAGS using:
#1. Linear spline regression
#1. Cubic B-spline
#2. Low-rank thin plate spline
#3. O'Sullivan spline


###########################################################
#1. Linear spline regression
X <- model.matrix(~ Lat.std + ML.std  + 
                 rhs(ML.std, -0.7121410) +  
                 rhs(ML.std, -0.1667513) +
                 rhs(ML.std, 0.6419299), 
                 data = Squid)
dim(X)
win.data <- list(Y = Squid$d15N,     
                 X = X,
                 M = ncol(X),          
                 N = nrow(Squid)
                  )
win.data

#####################################
#Model
sink("SquidGAM1.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M) { beta[i] ~ dnorm(0,0.0001) }
   
    #1B. Prior for variance for epsilon
    num   ~ dnorm(0, 0.0016) 
    denom ~ dnorm(0, 1)
    sigma <- abs(num / denom) 
    tau   <- 1 / (sigma * sigma)
       
    #2. Likelihood
    for (i in 1:N) {
       Y[i]  ~ dnorm(mu[i], tau)
       mu[i] <- inprod(X[i,], beta[])
              
       #3. Discrepancy measures   
       E[i] <- Y[i] - mu[i]
    }     
 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits  <- function () {
  list(beta  = rnorm(ncol(X), 0, 0.01),  #Regression parameters
       num   = rnorm(1, 0, 25),          #Prior stuff for variance epsilon
       denom = rnorm(1, 0, 1)            #Prior stuff for variance epsilon
       )  }


#Parameters to estimate
params <- c("beta", "E", "sigma")

#Start Gibbs sampler
K1   <- jags(data       = win.data,
             inits      = inits,
             parameters = params,
             model      = "SquidGAM1.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   = 14000,
             n.iter     = 15000)

K2 <- update(K1, n.iter = 20000) 
#Takes about 10 seconds on a new (2014) MacBook Pro.
print(K2, digits = 2)

out <- K2$BUGSoutput


#Not in the book...but an essential step
MyBUGSChains(out, c(uNames("beta", ncol(X)), "sigma"))
MyBUGSACF(out, uNames("beta", ncol(X)))
MyBUGSHist(out, uNames("beta", ncol(X)))

OUT1 <- MyBUGSOutput(out, c(uNames("beta", ncol(X)), "sigma"))
print(OUT1, digits =3)


#You will need to do some Googling to install this package
#Figure 1.6
library(coefplot2)
#lm
beta5 <- coef(M5)[2:6]
se5   <- sqrt(diag(vcov(M5)[2:6, 2:6]))

#JAGS 
beta1 <- OUT1[2:6,1]
se1   <- OUT1[2:6,2]



par(mar = c(5,5,2,2))
coefplot2(beta5, se5, offset = 0, col =1, 
          varnames = names(beta5), xlim = c(-2,3),
          cex.lab = 1.5, cex.var = 1)

coefplot2(beta1, se1, offset = 0.15, col = 1, 
          varnames = names(beta5), add = TRUE)
          
########################################################



#Access posterior mean residuals:
E <- out$mean$E


#######################################################
#Sketch smoother: Figure 1.7

#Extract coefficients
beta <- out$sims.list$beta  #ML smoother
dim(beta)

range(Squid$ML.std)
NewData <- data.frame(Lat.std = 0,
                      ML.std =  seq( -1.4,  3, length = 100))

Xnew <- model.matrix(~ Lat.std + ML.std  + 
                       rhs(ML.std, -0.7121410) +  
                       rhs(ML.std, -0.1667513) +
                       rhs(ML.std, 0.6419299), 
                     data = NewData)

#Calculate the smoothers
f <-  Xnew[,3:6] %*% t(beta[,3:6])


MySmoother <- function(f) {
 #Get the posterior mean, and 95% credible interval
 Sm1 <- matrix(nrow = 100, ncol = 5)
 for (i in 1:100) {
	Sm1[i,1:3] <- quantile(f[i,], 
                        probs = c(0.025, 0.5, 0.975))
    Sm1[i,4]  <- mean(f[i,])
	Sm1[i,5]  <- sd(f[i,])
	} 	
 colnames(Sm1) <- c("Lower", "Median", "Upper", "Mean", "SE")
 Sm1
 }
 
f.info <- MySmoother(f)


par(mar = c(5, 5, 2, 2))  
plot(x = NewData$ML.std, 
     y = f.info[,4], 
     type = "l", 
     xlab = "Standardized ML",
     ylab = "Smoother",
     ylim = c(-3, 2),
     cex.lab= 1.5, 
     lwd = 2)
lines(NewData$ML.std, f.info[,1], lty=2, lwd = 2)
lines(NewData$ML.std, f.info[,3], lty=2, lwd = 2)

################################################






########################################################################
#Section 1.5. B-splines in JAGS

I1      <- order(Squid$ML.std)
Squid2 <- Squid[I1,]
X      <- model.matrix(~Lat.std, data = Squid2)    

library(splines)
X.bs  <- bs(Squid2$ML.std, knots = probs[2:4], intercept = FALSE)
K     <- ncol(X.bs)

win.data2 <- list(Y            = Squid2$d15N,      #Response variable
                  X            = X,                #Covariates  
                  N            = nrow(Squid2),     #Sample size
                  X.bs         = X.bs,             #Basis for smoother 
                  M            = ncol(X),
                  K            = ncol(X.bs)
                  )
win.data2




#####################################
#Model
sink("SquidGAM2.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M) {beta[i] ~ dnorm(0, 0.0001)}  
    for (i in 1:K) {b[i]    ~ dnorm(0, 0.0001)}  
   
    #1B. Prior for variance for epsilon
    num   ~ dnorm(0, 0.0016) 
    denom ~ dnorm(0, 1)
    sigma <- abs(num / denom) 
    tau   <- 1 / (sigma * sigma)
       
    #2. Likelihood
    for (i in 1:N) {
       Y[i]  ~ dnorm(mu[i], tau)
       mu[i] <- inprod(X[i,], beta[]) + F1[i]
       F1[i] <- inprod(X.bs[i,], b[])  
              
       #3. Discrepancy measures   
       E[i] <- Y[i] - mu[i]
    }     
 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits2  <- function () {
  list(beta  = rnorm(ncol(X), 0, 0.01),        #Regression parameters
       b     = rnorm(ncol(X.bs), 0, 0.01), #Regression terms smoother
       num   = rnorm(1, 0, 25),          #Prior stuff for variance epsilon
       denom = rnorm(1, 0, 1)            #Prior stuff for variance epsilon
       )  }


#Parameters to estimate
params2 <- c("beta", "E", "sigma", "b", "F1")

#Start Gibbs sampler
K1   <- jags(data       = win.data2,
             inits      = inits2,
             parameters = params2,
             model      = "SquidGAM2.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   = 14000,
             n.iter     = 15000)

K2 <- update(K1, n.iter = 20000) 
print(K2, digits = 2)

out <- K2$BUGSoutput





MyBUGSChains(out, c(uNames("beta", ncol(X)), uNames("b", ncol(X.bs))))
#MyBUGSACF(out, uNames("beta", K))
#MyBUGSHist(out, uNames("beta", K))

OUT1 <- MyBUGSOutput(out, c(uNames("beta", ncol(X))))
print(OUT1, digits =3)

#Residual plots
E <- out$mean$E

################################################
#Sketch smoothers

#Extract coefficients
b <- out$sims.list$b
dim(b)

NewData <- data.frame(Lat.std = 0,
                      ML.std =  seq( -1.4,  3, length = 100))

X.bs  <- bs(NewData$ML.std, knots = probs[2:4], intercept = FALSE)

#Calculate the smoothers
f <-  X.bs %*% t(b)
#--->>>Why is the first row equal to 0????


MySmoother <- function(f) {
 #Get the posterior mean, and 95% credible interval
 Sm1 <- matrix(nrow = 100, ncol = 5)
 for (i in 1:100) {
	Sm1[i,1:3] <- quantile(f[i,], 
                        probs = c(0.025, 0.5, 0.975))
    Sm1[i,4]  <- mean(f[i,])
	Sm1[i,5]  <- sd(f[i,])
	} 	
 colnames(Sm1) <- c("Lower", "Median", "Upper", "Mean", "SE")
 Sm1
 }
 
f.info <- MySmoother(f)

 
 
#Not in the book 
#Plot the smoothers
par(mar = c(5,5,2,2))  
plot(x = NewData$ML.std, 
     y = f.info[,4], 
     type = "l", 
     xlab = "Standardized ML",
     ylab = "Smoother",
     ylim = c(-1,4),
     cex.lab= 1.5, lwd = 2)
lines(NewData$ML.std, f.info[,1], lty=2, lwd = 2)
lines(NewData$ML.std, f.info[,3], lty=2, lwd = 2)
#Odd!!! Better not do this?
##########################################################




##########################################################
#Section 1.6 Low rank thin plate regression spline in JAGS

K <- 5
Knots <- default.knots(Squid$ML.std, K)
Knots

#Figure 1.8

par(mar = c(5,5,2,2))
plot(y = Squid$d15N, 
     x = Squid$ML.std, 
     xlab = "Standardized ML",
     ylab = "d15N", cex.lab = 1.5,
     type = "n",
     ylim = c(0, 15))
abline(v = as.numeric(Knots), lty= 2)

x <- Squid$ML.std
x1 <- seq(from = min(x), max(x), length = 25)
for (i in 1:K){
  fi <- abs(x1-Knots[i])^3
 lines(x1,fi)
 }
  
text(1.25,12,"k = 1", cex = 1.2) 
text(1.72,12,"k = 2", cex = 1.2) 
text(2.75, 14,"k = 3", cex = 1.2) 
text(2.7, 10,"k = 4", cex = 1.2) 
text(2.7, 3,"k = 5", cex = 1.2) 

##############################################################




##############################################################

#Section 1.6.1
K <- 5
Knots <- default.knots(Squid$ML.std, K)
Z     <- GetZ_LRTP(Squid$ML.std, Knots)
X     <- model.matrix(~Squid$Lat.std + Squid$ML.std, data = Squid) 


dat1 <- data.frame(y = Squid$d15N, 
                   X = X, 
                   Z = Z, 
                   g = 1)
dat1$g <- factor(dat1$g)

library(nlme)
M1 <- lme(y~-1+X, random=list(g=pdIdent(~ Z-1)), data = dat1)
summary(M1)

beta <- M1$coef$fixed
u    <- unlist(M1$coef$random)
f    <- X[,3] * beta[3] + Z %*% u
varf <- X %*% vcov(M1) %*% t(X) 

#Smoother (not presented in the book)
I <- order(Squid$ML.std)
par(mar = c(5,5,2,2))
plot(x=sort(Squid$ML.std), f[I], type = "l", cex.lab = 1.5,
     xlab = "Standardized ML",
     ylab = "Smoother for ML")

sig.eps <- M1$sigma
sig.u   <- sig.eps * exp(unlist(M1$modelStruct))
####################################################################






####################################################################
#Section 1.6.2

#Fixed covariates
Xcov <- model.matrix(~Lat.std, data = Squid)    

#Get the X and Z for the tprs
#Knots:
K <- 5
Knots <- default.knots(Squid$ML.std, K)
Z <- GetZ_LRTP(Squid$ML.std, Knots)

N <- nrow(Squid)
X <- Squid$ML.std

win.data2 <- list(Y            = Squid$d15N,          #Response variable
                  Xcov         = Xcov,                #Covariates  
                  N            = nrow(Squid),         #Sample size
                  M            = ncol(Xcov),          #Number of betas
                  Xspl         = X,                   #Basis for smoother 
                  Zspl         = Z,
                  Mz           = ncol(Z)              #Number of us
                  )
win.data2


#We are using multivariate Normal priors here..and also in the book.
#Having see a lot of other examples during our GAMM courses (where we have
#20 participants all running this, and other data sets), I think using
#univariate Normal priors works better.

#####################################
#Model
sink("GAMlrtprs.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M)  {beta[i] ~ dnorm(0, 0.0001) }  
    for (i in 1:Mz) {u[i] ~ dnorm(0, tau.u )  }
    b ~ dnorm(0, 0.0001)
   
    #1B. Prior for variance for epsilon
    num   ~ dnorm(0, 0.0016) 
    denom ~ dnorm(0, 1)
    sigma <- abs(num / denom) 
    tau   <- 1 / (sigma * sigma)
       
    #Priors for variance random intercept u for smoother
    num.u   ~ dnorm(0, 0.0016) 
    denom.u ~ dnorm(0, 1)
    sigma.u <- abs(num.u / denom.u) 
    tau.u   <- 1 / (sigma.u * sigma.u)
     
    #2. Likelihood
    for (i in 1:N) {
       Y[i]  ~  dnorm(mu[i], tau)
       mu[i] <- inprod(beta[], Xcov[i,]) + F1[i]         
       F1[i] <- Xspl[i] * b + inprod(Zspl[i,], u[])  #Smoother
              
       #3. Discrepancy measures   
       E[i] <- Y[i] - mu[i]
    }     
 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits2  <- function () {
  list(beta    = rnorm(ncol(Xcov), 0, 0.1),#Regression parameters
       b       = rnorm(1, 0, 0.1),         #Regression terms smoother
       num     = rnorm(1, 0, 25),          #Prior stuff for variance epsilon
       denom   = rnorm(1, 0, 1),           #Prior stuff for variance epsilon
       num.u   = rnorm(1, 0, 25),          #Prior stuff for variance u
       denom.u = rnorm(1, 0, 1),           #Prior stuff for variance u
       u       = rnorm(K, 0, 1)
       )  }


#Parameters to estimate
params2 <- c("beta", "E", "sigma", "b", "u", "sigma.u")

#Start Gibbs sampler
K1   <- jags(data       = win.data2,
             inits      = inits2,
             parameters = params2,
             model      = "GAMlrtprs.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   =  15000,
             n.iter     = 115000)

#K2 <- update(K1, n.iter = 100000) 

print(K1, digits = 2)
out <- K1$BUGSoutput



MyBUGSChains(out, c(uNames("beta", ncol(Xcov))))
MyBUGSChains(out, c("b",uNames("u", ncol(Z))))
MyBUGSChains(out, c("sigma", "sigma.u"))

##MyBUGSACF(out, uNames("beta", K))
##MyBUGSHist(out, uNames("beta", K))

OUT1 <- MyBUGSOutput(out, c(uNames("beta", ncol(Xcov))))
print(OUT1, digits =3)


################################################
#Sketch smoothers

#Extract coefficients
b <- out$sims.list$b  #ML smoother
u <- out$sims.list$u


#This depends on the selected data!
range(Squid$ML.std)
ML.100 <- seq( -1.4,  3, length = 100)

#Convert this covariate into a smoother basis

Z100 <- GetZ_LRTP(ML.100, Knots)


#Calculate the smoothers
f1 <-  ML.100 %*% t(b) + Z100 %*% t(u) 

#Get posterior mean and 95% credible intervals
f1.info <- MySmoother(f1)

 
#Figure 1.9 

#Plot the smoothers
par(mar = c(5,5,2,2))  
plot(x = ML.100, 
     y = f1.info[,4], 
     type = "l", 
     xlab = "Standardized ML",
     ylab = "Smoother",
     ylim = c(-2,2), cex.lab = 1.5)
lines(ML.100, f1.info[,1], lty=2)
lines(ML.100, f1.info[,3], lty=2)


#######################################################################
#Figuer 1.10
#This will take a while!!!
#Write a loop in which we use K = 2, 3,4,5,6,7,8,9,10,11,12,13,14,15
fm  <- NULL
fl  <- NULL
fu  <- NULL
fid <- NULL

for (k in 2:15) {
  K <- k
  Knots <- default.knots(Squid$ML.std, K)
  Z <- GetZ_LRTP(Squid$ML.std, Knots)

  win.data2 <- list(Y    = Squid$d15N,     #Response variable
                    Xcov = Xcov,           #Covariates  
                    N    = nrow(Squid),    #Sample size
                    M    = ncol(Xcov),     #Betas 
                    Xspl = X,              #Basis for smoother 
                    Zspl = Z,
                    Mz   = ncol(Z)         #Number of us for smoother basis
                  )

 K1   <- jags(data       = win.data2,
              inits      = inits2,
              parameters = params2,
              model      = "GAMlrtprs.txt",
              n.thin     = 10,
              n.chains   = 3,
              n.burnin   = 14000,
              n.iter     = 15000)

  K2 <- update(K1, n.iter = 20000) 
  out <- K2$BUGSoutput

  b <- out$sims.list$b  
  u <- out$sims.list$u
  ML.100 <- seq( -1.4,  3, length = 100)
  Z100 <- GetZ_LRTP(ML.100, Knots)
  f1 <-  ML.100 %*% t(b) + Z100 %*% t(u) 
  f1.info <- MySmoother(f1)

  fm <- c(fm,f1.info[,4])
  fl <- c(fl,f1.info[,1])
  fu <- c(fu,f1.info[,3])
  fid<- c(fid,rep(k, 100))
  print(k)
}

fx <- rep(ML.100, 14)


xyplot(fm ~ fx | factor(fid),
    xlab = list(label = "Standardized ML", cex = 1.5),
    ylab = list(label = "Smoother", cex = 1.5),
    ylim = c(-2,2),
    panel = function(x, y, subscripts) {
    	panel.lines(x, y, lwd = 2, col = 1)
    	panel.lines(x, fl[subscripts], lty = 2, lwd = 2, col = 1)
    	panel.lines(x, fu[subscripts], lty = 2, lwd = 2, col = 1)
    })
    
####################################################################






######################################################
#Section 1.7. O'Sullivan splines in JAGS

Xcov <- model.matrix(~Lat.std, data = Squid)    
numIntKnots <- 11
probs <- seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))]
intKnotsTime <- quantile(unique(Squid$ML.std), probs)

XZML <- OSullivan(Squid$ML.std, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsTime)

#As mentioned earlier in this file, I think that using
#univariate Normal priors may be better then multivariate Normal
#priors.
win.data3 <- list(Y      = Squid$d15N,    #Response variable
                  Xcov   = Xcov,          #Covariates  
                  N      = nrow(Squid),   #Sample size
                  M      = ncol(Xcov),    #Betas 
                  Xspl   = XZML$X[,1],    #Basis for smoother 
                  Zspl   = XZML$Z,
                  Mz     = ncol(XZML$Z)   #u for smoother basis
                 )
win.data3

#Model
sink("GAMOSS.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M)  { beta[i] ~ dnorm(0, 0.0001) }  
    for (i in 1:Mz) { u[i] ~ dnorm(0, tau.u) }  
    b ~ dnorm(0, 0.0001)
   
    #1B. Prior for variance for epsilon
    num   ~ dnorm(0, 0.0016) 
    denom ~ dnorm(0, 1)
    sigma <- abs(num / denom) 
    tau   <- 1 / (sigma * sigma)
       
    #Priors for variance random intercept u for smoother
    num.u   ~ dnorm(0, 0.0016) 
    denom.u ~ dnorm(0, 1)
    sigma.u <- abs(num.u / denom.u) 
    tau.u   <- 1 / (sigma.u * sigma.u)
        
    #2. Likelihood
    for (i in 1:N) {
       Y[i]   ~ dnorm(mu[i], tau)
       mu[i] <- inprod(beta[], Xcov[i,]) + F1[i]         
       F1[i] <- Xspl[i] * b + inprod(Zspl[i,], u[]) 
              
       #3. Discrepancy measures   
       E[i] <- Y[i] - mu[i]
    }     
 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits3  <- function () {
  list(beta  = rnorm(ncol(Xcov), 0, 0.1),  #Regression parameters
       b     = rnorm(1, 0, 0.01),
       u     = rnorm(ncol(XZML$Z), 0, 0.1),#Regression terms smoother
       num   = rnorm(1, 0, 25),            #Prior stuff for variance epsilon
       denom = rnorm(1, 0, 1),             #Prior stuff for variance epsilon
       num.u   = rnorm(1, 0, 25),          #Prior stuff for variance u
       denom.u = rnorm(1, 0, 1)            #Prior stuff for variance u
       )  }


#Parameters to estimate
params3 <- c("beta", "E", "sigma", "b", "u", "sigma.u")


#Start JAGS. Takes about 1 minute in a 2014 Macbook Pro
K1   <- jags(data       = win.data3,
             inits      = inits3,
             parameters = params3,
             model      = "GAMOSS.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   = 15000,
             n.iter     = 115000)

out <- K1$BUGSoutput


#Assess mixing
MyBUGSChains(out, c(uNames("beta", ncol(Xcov)), uNames("u", 3), "b", "sigma", "sigma.u"))
MyBUGSACF(out, uNames("beta", ncol(Xcov)))
MyBUGSHist(out, uNames("beta", ncol(Xcov)))

OUT1 <- MyBUGSOutput(out, c(uNames("beta", ncol(Xcov)),"sigma", "sigma.u"))
print(OUT1, digits =3)

#Use for residual plots:
E <- out$mean$E

################################################
#Figure 1.11
#Extract coefficients
b <- out$sims.list$b  
u <- out$sims.list$u


#This depends on the selected data!
range(Squid$ML.std)
ML.100 <- seq( -1.4,  3, length = 100)

#Convert this covariate into a smoother basis

XZML100 <- OSullivan(ML.100, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsTime)

Z100 <- XZML100$Z
#Calculate the smoothers
f1 <-  ML.100 %*% t(b) + Z100 %*% t(u) 

#Get posterior mean and 95% credible intervals
f1.info <- MySmoother(f1)

 
#Plot the smoothers



par(mar = c(5,5,2,2), cex.lab = 1.5)  
plot(x = ML.100, 
     y = f1.info[,4], 
     type = "l", 
     xlab = "Standardized ML",
     ylab = "Smoother",
     ylim = c(-2,2), cex.lab = 1.5)
lines(ML.100, f1.info[,1], lty=2)
lines(ML.100, f1.info[,3], lty=2)



##################################
#Section 1.8 Degrees of freedom
XFixed <- model.matrix(~ 1 + Lat.std, data = Squid)    

numIntKnots <- 11
probs <- seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))]
intKnotsTime <- quantile(unique(Squid$ML.std), probs)

XZML <- OSullivan(Squid$ML.std, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsTime)

C <- cbind(XFixed, win.data3$Xspl, win.data3$Z)

D <- diag(ncol(win.data3$Z)+3)
D[1,1] <- 0
D[2,2] <- 0
D[3,3] <- 0


lambda <- out$mean$sigma^2 / out$mean$sigma.u^2
lambda <- as.numeric(lambda) 
bu     <- solve(t(C) %*% C + lambda * D) %*% t(C) %*% win.data3$Y
Fit    <- C %*% bu





par(mar = c(5,5,2,2), cex.lab = 1.5)  
plot(x = ML.100, 
     y = f1.info[,4], 
     type = "l", 
     xlab = "Standardized ML",
     ylab = "Smoother",
     ylim = c(-2,2), cex.lab = 1.5)
lines(ML.100, f1.info[,1], lty=2)
lines(ML.100, f1.info[,3], lty=2)


#Run this code after the plot command for the smoother in Equation (1.20)
f_lambda <- cbind(XZML100$X, XZML100$Z) %*% bu[3:16]
lines(ML.100, f_lambda, lty=2, col = 1, lwd = 2)
legend("topleft",
       legend = c("Smoother based on MCMC",
                  "Smoother based on Equation (1.20)"),
       lwd = c(1,2),
       lty = c(1,2),
       cex = 1.2)


#Degrees of freedom
Something <- C %*% solve(t(C) %*% C + lambda * D) %*% t(C)
sum(diag(Something))

Someth <- solve(t(C) %*% C + lambda * D) %*% t(C)
xx <- C[,3:16] %*% Someth[3:16,]
sum(diag(xx))



