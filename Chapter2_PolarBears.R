# R code for Chapter 2 in: Beginner's Guide to GAMM with R (2014)
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
setwd("/Users/Highstat/applicat/HighlandStatistics/Books/BGS/GAMM/Data/PolarBears")


#Load the data
PB <- read.table(file="PolarBearsV2.txt", header=TRUE)

#Check what we have:
names(PB)
str(PB)
#############################################################################



#####################################################################
#Load packages
library(lattice)
library(mgcv)
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV6.R")  
source("/Users/Highstat/applicat/HighlandStatistics/MCMC/R/SupportFilesHighStat.R")
source(file = "/Users/Highstat/applicat/HighlandStatistics/MCMC/R/MCMCSupportHighstat.R")
library(R2jags)
##################################################################





######################################################################
#House keeping
PB$fBearID   <- factor(PB$BearID)
PB$fRepro    <- factor(PB$Repro)
PB$fDen      <- factor(PB$Den)
PB$fMonth    <- factor(PB$Month)
PB$fSeason   <- factor(PB$Season)
######################################################################




######################################################################
#Data exploration
#Outliers
dotchart(PB$Movement, 
         xlab = "Values of the data",
         ylab = "Order of the data")

#Graph for Word (takes less memory)
#Figure 2.2
par(mar=c(5,5,2,2), cex.lab = 1.5)
plot(y = 1:nrow(PB), 
     x = PB$Movement, 
     xlab = "Values of the data",
     ylab = "Order of the data",
     pch = 16, 
     cex = 0.7)


#Relationships
#Figure 2.3
par(mar=c(5,5,2,2), cex.lab = 1.5)
plot(y = PB$Movement, 
     x = PB$dDay,
     xlab = "Days since 1 January 1988",
     ylab = "Movement")

#Add a smoother          
M1      <- loess(Movement ~ dDay, data = PB)
MyData1 <- data.frame(dDay = seq(98,4362))
P1      <- predict(M1, newdata = MyData1)
lines(x = MyData1$dDay, 
      y = P1, 
      lwd = 5)



#Figure 2.4
par(mar=c(5,5,2,2), cex.lab = 1.5)
boxplot(Movement ~ fMonth,
        data = PB,
        xlab = "Month", 
        ylab = "Movement",
        varwidth = TRUE)
#Seems to be more travelling during the winter

#Figure 2.5
xyplot(Movement ~ DayInYear | factor(Year),
       col = 1, 
       data = PB,
       cex=0.5,
       xlab = list("Day in year", cex = 1.5),
       ylab = list("Movement", cex = 1.5),
       strip = function(bg='white',...) strip.default(bg='white', ...),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x,y,col=1)
       })


#Figure 2.6
xyplot(Movement ~ dDay | factor(BearID),
  col = 1, data = PB,
  cex=0.5,
  strip = function(bg='white',...) strip.default(bg='white', ...),
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")),
  xlab = "Day number",
  ylab = "Movement",
  panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    })



#How many bears where tagged simultanously?
#Figure 2.7
MF <- function(x){ length(unique(x)) }

NBears <- tapply(PB$BearID, 
                 INDEX = PB$Year, 
                 FUN = MF)

par(mar=c(5,5,2,2), cex.lab = 1.5)
plot(x = 1988:1999, 
     y = NBears, 
     pch = 16,
     xlab = "Year",
     ylab = "Number of bears with working tag")



#Figure 2.8
par(mar=c(5,5,2,2), cex.lab = 1.5)
boxplot(Movement ~ BearID, 
        data = PB,
        xlab = "Polar bear",
        ylab = "Movement",
        varwidth = TRUE)



#Figure 2.9
#Do we have data from the same months over time?
par(mar=c(5,5,2,2), cex.lab = 1.5)
plot(x = jitter(PB$Year), 
     y = jitter(as.numeric(PB$fMonth)),
     xlab = "Year",
     ylab = "Month")
     
M1      <- loess(as.numeric(fMonth) ~ Year, data = PB)
MyData1 <- data.frame(Year = seq(1988,1999))
P1      <- predict(M1, newdata = MyData1)
lines(x = MyData1$Year, 
      y = P1, 
      lwd = 5)



table(PB$fDen, PB$fRepro)
#Not good for an interaction!
###################################################






#Section 2.6
PB$Yearc      <- PB$Year - mean(PB$Year)
PB$DayInYearc <- PB$DayInYear - mean(PB$DayInYear)

M1 <- gamm(Movement ~ fRepro + Yearc + s(DayInYearc),
           random = list(BearID=~1),
           data = PB, method = "ML")

summary(M1$gam)
anova(M1$gam)
summary(M1$lme)


#Figure 2.10
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(M1$gam, cex.lab = 1.5, axes = FALSE, xlab = "Day in year")
axis(2)


AxisScaled <- seq(from = min(PB$DayInYearc),
                  to = max(PB$DayInYearc), 
                  length = 5)
AxisOriginal <- AxisScaled + mean(PB$DayInYear)
axis(1, at = AxisScaled, labels = AxisOriginal)
box()

        
#Figure 2.11    
par(mfrow = c(2,2), mar = c(5,5,2,2), cex.lab = 1.5)    
E1 <- resid(M1$lme, type ="n")
F1 <- fitted(M1$lme)
plot(x=F1, y=E1, xlab = "Fitted values", ylab ="Residuals")
abline(h=0, lty=2)

plot(x=PB$Year, y = E1, xlab = "Year", ylab = "Residuals")
abline(h=0, lty=2)

plot(x=PB$DayInYear, y = E1, xlab = "Day in Year", ylab = "Residuals")
abline(h=0, lty=2)

boxplot(E1 ~ Season, data = PB)



M2 <- gamm(Movement ~ fRepro + Yearc + s(DayInYearc),
           random = list(BearID=~1),
           weight = varIdent(form = ~1 | fSeason),
           data = PB, REML = TRUE)

summary(M2$gam)
summary(M2$lme)
anova(M2$gam)
plot(M2$gam)


M1 <- gamm(Movement ~ fRepro + Yearc + s(DayInYearc),
           random = list(BearID=~1),
           data = PB, method = "REML")

M2 <- gamm(Movement ~ fRepro + Yearc + s(DayInYearc),
           random = list(BearID=~1),
           weight = varIdent(form = ~1 | fSeason),
           data = PB, method = "REML")

anova(M1$lme, M2$lme)


sigmaA <- 14.76934 
sigmaEps <- 33.54436
Sj <- c(1.000000, 1.518950, 1.132879, 1.854835)
sigmaA^2 / (sigmaA^2 + (Sj * sigmaEps)^2)
###########################################################################






###########################################################################
#Section 2.6.3

Xcov <- model.matrix(~1 + fRepro + Yearc, data = PB)

numIntKnots <- 5
probs <- seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))]
intKnotsDiY <- quantile(unique(PB$DayInYearc), probs)

XZDayInYear <- OSullivan(PB$DayInYearc, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsDiY)

X <- cbind(Xcov, XZDayInYear$X[,1])
Z <- XZDayInYear$Z
dat1 <- data.frame(y = PB$Movement, 
                   X = X, 
                   Z = Z, 
                   g = 1)
                                     
dat1$g <- factor(dat1$g)

library(nlme)
M3 <- lme(y~-1+X, 
         random=list(g=pdIdent(~ Z-1)), 
         data = dat1)

summary(M1$lme)         
summary(M3)


beta.hat <- M3$coef$fixed
u.hat <- unlist(M3$coef$random)
f.hat <- X[,5] * beta.hat[5] + Z %*% u.hat

sig.eps.hat <- M3$sigma
sig.u.hat <- sig.eps.hat*exp(unlist(M3$modelStruct))

#Figure 2.12
I <- order(PB$DayInYear)
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(x=sort(PB$DayInYear), f.hat[I], type = "l", cex.lab = 1.5,
     xlab = "Day in Year",
     ylab = "Smoother for Day in Year")
######################################################################







######################################################################
#Section 2.7 
Xcov <- model.matrix(~1 + Yearc + fRepro , data = PB)

#Random effect stuff
re    <- as.numeric(as.factor(PB$BearID))
NumPB <- length(unique(PB$BearID))


#Spline
numIntKnots <- 5
probs <- seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))]
intKnotsDiY <- quantile(unique(PB$DayInYearc), probs)

XZDayInYear <- OSullivan(PB$DayInYearc, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsDiY)


#Allow for 4 variances
as.numeric(PB$fSeason)

#As above...better use univariate Normal priors
win.data1 <- list(Y            = PB$Movement,            #Response variable
                  Xcov         = Xcov,                   #Covariates  
                  N            = nrow(PB),               #Sample size
                  M            = ncol(Xcov),             #Number of betas 
                  Xspl         = XZDayInYear$X[,1],      #Basis for smoother 
                  Zspl         = XZDayInYear$Z,
                  Mz           = ncol(XZDayInYear$Z),    #Number of us for smoother basis
                  re           = re,                     #Random effects
                  NumPB        = NumPB,                  #Number of random effects
                  Season       = as.numeric(PB$fSeason)
                  )
win.data1



#####################################
#Model
sink("AMOSS.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M) {beta[i] ~ dnorm(0, 0.0001) }  
      
    #1B. Prior for variance for epsilon
    for (i in 1:4) {
      num[i]   ~ dnorm(0, 0.0016) 
      denom[i] ~ dnorm(0, 1)
      sigma[i] <- abs(num[i] / denom[i]) 
      tau[i]   <- 1 / (sigma[i] * sigma[i])
    }
    
    #1C. Priors for smoother components
    for (i in 1:Mz) { u[i] ~ dnorm(0, tau.u) } 
    b ~ dnorm(0, 0.0001)

    #1D. Priors for variance random intercept u for smoother
    num.u   ~ dnorm(0, 0.0016) 
    denom.u ~ dnorm(0, 1)
    sigma.u <- abs(num.u / denom.u) 
    tau.u   <- 1 / (sigma.u * sigma.u)
     
    #1E. Priors for random effects polar bear
    for (i in 1:NumPB) {a[i] ~ dnorm(0, tau.re) }
    
    #1F. Priors for sigma for random effects polar bear
    num.re   ~ dnorm(0, 0.0016) 
    denom.re ~ dnorm(0, 1)
    sigma.re <- abs(num.re / denom.re) 
    tau.re   <- 1 / (sigma.re * sigma.re)
 
    #2. Likelihood
    for (i in 1:N) {
       Y[i]  ~ dnorm(mu[i], tau[Season[i]])
       mu[i] <- inprod(beta[], Xcov[i,]) + F1[i]  + a[re[i]]       
       F1[i] <- Xspl[i] * b + inprod(Zspl[i,], u[]) 
              
       #3. Discrepancy measures   
        YNew[i]   ~ dnorm(mu[i], tau[Season[i]]) 
        expY[i]    <- mu[i] 
        varY[i]    <- sigma[Season[i]] * sigma[Season[i]]
        E[i]       <- Y[i] - mu[i]
        PRes[i]    <- (Y[i]  - expY[i]) / sqrt(varY[i])
        PResNew[i] <- (YNew[i] - expY[i]) / sqrt(varY[i])
        D[i]       <- pow(PRes[i], 2)
        DNew[i]    <- pow(PResNew[i], 2)
    }     
     Fit         <- sum(D[1:N])
     FitNew      <- sum(DNew[1:N])        
 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits1  <- function () {
  list(beta      = rnorm(ncol(Xcov), 0, 0.1),        #Regression parameters
       a         = rnorm(NumPB, 0, 0.1),
       num.re    = rnorm(1, 0, 25),   
       denom.re  = rnorm(1, 0, 1),
       b         = rnorm(1, 0, 0.1),
       u         = rnorm(ncol(XZDayInYear$Z), 0, 0.01), #Regression terms smoother
       num   = rnorm(4, 0, 25),          #Prior stuff for variance epsilon
       denom = rnorm(4, 0, 1),           #Prior stuff for variance epsilon
       num.u   = rnorm(1, 0, 25),        #Prior stuff for variance u
       denom.u = rnorm(1, 0, 1)          #Prior stuff for variance u
       )  
}

#Parameters to estimate
params1 <- c("beta", "E", "sigma", "b", "u", 
             "sigma.u", "sigma.re", "Fit", "FitNew" )

#Start Gibbs sampler
# K1   <- jags(data       = win.data1,
             # inits      = inits1,
             # parameters = params1,
             # model      = "AMOSS.txt",
             # n.thin     = 10,
             # n.chains   = 3,
             # n.burnin   = 15000,
             # n.iter     = 25000)
#out <- K1$BUGSoutput


#You better use this
K1   <- jags(data       = win.data1,
             inits      = inits1,
             parameters = params1,
             model      = "AMOSS.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   = 4000,
             n.iter     = 5000)
K2  <- update(K1, n.iter = 10000, n.thin = 10)  
out <- K2$BUGSoutput
out

#Check mixing
#Figure 2.13
MyBUGSChains(out, c(uNames("beta", ncol(Xcov)), 
                    uNames("u", ncol(XZDayInYear$Z)), 
                    "b", 
                    "sigma.u", "sigma.re",
                    uNames("sigma", 4)))

#Not in the book:
MyBUGSACF(out, uNames("beta", ncol(Xcov)))
MyBUGSHist(out, uNames("beta", ncol(Xcov)))



#Section 2.7.7
E      <- out$mean$E
Sigma  <- out$mean$sigma[PB$Season]
E.norm <- E / Sigma


OUT1 <- MyBUGSOutput(out,  c(uNames("beta", ncol(Xcov)), 
                    uNames("u", ncol(XZDayInYear$Z)), 
                    "b", 
                    "sigma.u", "sigma.re",
                    uNames("sigma", 4)))
                    
print(OUT1, digits = 2)


#Sketch smoother

#Extract coefficients
b <- out$sims.list$b  
u <- out$sims.list$u


#This depends on the selected data!
range(PB$DayInYearc)
DayInYearc.100 <- seq( -167.43,  196.56, length = 100)

#Convert this covariate into a smoother basis

XZ100 <- OSullivan(DayInYearc.100, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsDiY)

Z100 <- XZ100$Z
#Calculate the smoothers
f1 <-  XZ100$X %*% t(b) + Z100 %*% t(u) 
dim(f1)
#Get posterior mean and 95% credible intervals
f1.info <- MySmoother(f1)

 
#Plot the smoothers

STDScale <- quantile(DayInYearc.100, c(0, 0.5, 1))
OriScale <- STDScale * 1 + mean(PB$DayInYear)
cbind(OriScale, STDScale)

DayInYear.100 <- DayInYearc.100 * 1 + mean(PB$DayInYear)

par(mar = c(5,5,2,2), cex.lab = 1.5)  
plot(x = DayInYear.100, 
     y = f1.info[,4], 
     type = "l", 
     xlab = "Day in Year",
     ylab = "Smoother",
     ylim = c(-25,25), cex.lab = 1.5,
     axes = TRUE)
#axis(2)
#axis(1, at = STDScale, labels = OrScale)
     
lines(DayInYear.100, f1.info[,1], lty=2)
lines(DayInYear.100, f1.info[,3], lty=2)
quartz()
plot(M2$gam)




############################################################################
#This is part of the model validation (and is not presented in the book)
par(mfrow = c(2,2), mar = c(5,5,2,2))    
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values", 
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = PB$Year, 
     y = E.norm, 
     xlab = "Year", 
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = PB$DayInYear, 
     y = E.norm, 
     xlab = "Day in Year", 
     ylab = "Residuals")
abline(h=0, lty=2)

boxplot(E.norm ~ Season, data = PB)
############################################################################




#############################################
#Gamma GAM
M3 <- gamm(Movement ~ fRepro + Yearc + s(DayInYearc),
           random = list(BearID=~1),
           family = Gamma(link = "log"),
           data = PB)

#Not in the book:
summary(M3$gam)
summary(M3$lme)
anova(M3$gam)
plot(M3$gam)

E3 <- resid(M3$lme, type = "pearson")
F3 <- fitted(M3$gam)
plot(x = F3, 
     y = E3)

E3a <- resid(M3$lme, type = "pearson")
E3b <- resid(M3$gam, type = "pearson")
#One type contains the random effects.

mu <- fitted(M3$lme)
EP <- (PB$Movement - mu) / sqrt(mu)

EP - E3b

a <- ranef(M3$lme)$BearID$'(Intercept)'
exp(fitted(M3$lme) + a[re])
fitted(M3$gam)
##################################################





##################################################
#JAGS Gamma GAMM
##################################################


win.data2 <- list(Y            = PB$Movement,         #Response variable
                  Xcov         = Xcov,                #Covariates  
                  N            = nrow(PB),            #Sample size
                  M            = ncol(Xcov),          #Number of betas 
                  Xspl         = XZDayInYear$X[,1],   #Basis for smoother 
                  Zspl         = XZDayInYear$Z,
                  Mz           = ncol(XZDayInYear$Z), #Number of us for smoother basis
                  re           = re,                  #Random effects
                  NumPB        = NumPB                #Number of RE bear
                  )
win.data2



#####################################
#Model (better use univariate Normal priors)
sink("GAMMOSS.txt")
cat("
model{
    #1A. Priors regression parameters
    for (i in 1:M) { beta[i] ~ dnorm(0, 0.001) }
      
    #1B. Prior for r parameter of Gamma distribution
    r ~ dgamma( 0.01, 0.01 )
    
    #1C. Priors for smoother components
    for (i in 1:Mz) { u[i] ~ dnorm(0, tau.u) }  
    b ~ dnorm(0, 0.0001)

    #1D. Priors for variance random intercept u for smoother
    num.u   ~ dnorm(0, 0.0016) 
    denom.u ~ dnorm(0, 1)
    sigma.u <- abs(num.u / denom.u) 
    tau.u   <- 1 / (sigma.u * sigma.u)
     
    #1E. Priors for random effects polar bear
    for (i in 1:NumPB) { a[i] ~ dnorm(0, tau.re) } 
    
    #1F. Priors for sigma for random effects polar bear
    num.re   ~ dnorm(0, 0.0016) 
    denom.re ~ dnorm(0, 1)
    sigma.re <- abs(num.re / denom.re) 
    tau.re   <- 1 / (sigma.re * sigma.re)
 
    #2. Likelihood
    for (i in 1:N) {
       Y[i]  ~ dgamma(r, mu.eff[i])
       mu.eff[i]  <- r / mu[i]
       log(mu[i]) <- eta[i]
       #mu[i] <- eta[i]
       eta[i] <- inprod(beta[], Xcov[i,]) + F1[i]  + a[re[i]]       
       F1[i] <- Xspl[i] * b + inprod(Zspl[i,], u[]) 
              
       #3. Discrepancy measures   
        YNew[i]   ~ dgamma(r, mu.eff[i]) 
        expY[i]    <- mu[i] 
        varY[i]    <- (1 / r) * mu[i]* mu[i]
        PRes[i]    <- (Y[i]  - expY[i]) / sqrt(varY[i])
        PResNew[i] <- (YNew[i] - expY[i]) / sqrt(varY[i])
        D[i]       <- pow(PRes[i], 2)
        DNew[i]    <- pow(PResNew[i], 2)
    }     
     Fit         <- sum(D[1:N])
     FitNew      <- sum(DNew[1:N])

 }
",fill = TRUE)
sink()
#####################################
#

 
#Inits function
inits2  <- function () {
  list(beta      = rnorm(ncol(Xcov), 0, 0.1),        #Regression parameters
       a         = rnorm(NumPB, 0, 0.1),
       num.re    = rnorm(1, 0, 25),   
       denom.re  = rnorm(1, 0, 1),
       b         = rnorm(1, 0, 0.01),
       u         = rnorm(ncol(XZDayInYear$Z), 0, 0.1), #Regression terms smoother
       r         = 1,
       num.u   = rnorm(1, 0, 25),          #Prior stuff for variance u
       denom.u = rnorm(1, 0, 1)            #Prior stuff for variance u
       )  }


#Parameters to estimate
params2 <- c("beta", "r", "b", "u", 
             "sigma.u", "sigma.re", "Fit", "FitNew", "mu", "PRes" )

#Start Gibbs sampler
K1   <- jags(data       = win.data2,
             inits      = inits2,
             parameters = params2,
             model      = "GAMMOSS.txt",
             n.thin     = 10,
             n.chains   = 3,
             n.burnin   = 4000,
             n.iter     = 5000)
#I prefer this:
K2   <- update(K1, n.iter = 10000, n.thin = 10)
out2 <- K2$BUGSoutput
out2

mean(K1$BUGSoutput$sims.list$Fit >  
     K1$BUGSoutput$sims.list$FitNew) 



MyBUGSChains(out, c(uNames("beta", ncol(Xcov)), 
                    uNames("u", 3), 
                    "r",
                    "b", "r",
                    "sigma.u", "sigma.re"))
#MyBUGSACF(out, uNames("beta", K))
#MyBUGSHist(out, uNames("beta", K))

OUT2 <- MyBUGSOutput(out2, c(uNames("beta", ncol(Xcov)),
                            "r"))
print(OUT2, digits =3)



PRes <- out2$mean$PRes
Fit  <- out2$mean$mu

plot(x = Fit, y = PRes)

#Sketch smoother

#Extract coefficients
b <- out2$sims.list$b  
u <- out2$sims.list$u

DayInYearc.100 <- seq( -167.43,  196.56, length = 100)
#Convert this covariate into a smoother basis

XZ100 <- OSullivan(DayInYearc.100, 
                    numIntKnots = numIntKnots, 
                    AddIntercept = FALSE,
                    intKnots = intKnotsDiY)

Z100 <- XZ100$Z
#Calculate the smoothers
f1 <-  XZ100$X %*% t(b) + Z100 %*% t(u) 

#Get posterior mean and 95% credible intervals
f1.info <- MySmoother(f1)

 
#Plot the smoothers
DayInYear.100 <- DayInYearc.100 * 1 + mean(PB$DayInYear)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)  
plot(x = DayInYear.100, 
     y = f1.info[,4], 
     type = "l", 
     xlab = "Day in Year",
     ylab = "Smoother",
     ylim = c(-0.4,0.5), cex.lab = 1.5,
     axes = TRUE)
     
lines(DayInYear.100, f1.info[,1], lty=2)
lines(DayInYear.100, f1.info[,3], lty=2)

plot(M3$gam, cex.lab = 1.5, axes = FALSE, xlab = "Day in year")
axis(2)
AxisScaled <- seq(from = min(PB$DayInYearc),
                  to = max(PB$DayInYearc), 
                  length = 5)
AxisOriginal <- AxisScaled + mean(PB$DayInYear)
axis(1, at = AxisScaled, labels = AxisOriginal)
box()
