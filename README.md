# Efficacy-of-Location-level-Fixed-Effects-in-Experimental-Economics

>Note: This repository corresponds to the following pre-print: (TO FILL: Insert arXiv link)

### What the R Code Does
- Simulates data with heterogeneity in **location**
- Tests nine meta-regression methodologies
- Calculates estimator bias, MSE, power, and CI coverage
- Generates figures summarizing results

### Key Workflow
- Vary number of countries: 5, 19, 33
- Vary location heterogeneity
- Repeat simulations via Monte Carlo design

### Metrics Output
- Power, RMSE, bias, confidence intervals, coverage rates
- Exported as data frames for comparison

### Sample Code

This illustrates how the code operates when there is one covariate

#### Set Up - Find True Parameters
```

### Start by setting up variance covariance matrix (presumed same for all studies)
library(mvtnorm)
library(Matrix)

# matrix syntax: (TL, BL, TR, BR)
S00 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)      # top left
S01 <- matrix(c(0.25, 0, 0.2, 0.1), nrow = 2, ncol = 2)   # top right
S10 <- t(S01)                                             # bottom left
S11 <- matrix(c(1, 0.1, 0.1, 1), nrow = 2, ncol = 2)      # bottom right

SIGMA <- rbind(cbind(S00, S01), cbind(S10, S11))          # VarCov matrix

# Elements used to calculate true parameters
varY <- S00[1,1]
covY <- SIGMA[2:length(SIGMA[1,]),1:1]
covX <- SIGMA[2:length(SIGMA[1,]),2:length(SIGMA[1,])]   
# the length() extracts the length of one row
trueBetas <- solve(covX) %*% t(t(covY)) # t(t()) to transpose once, likely because of datatype

### The true parameters
sig2True <- varY - covY %*% trueBetas                     # Omega
b1True <- trueBetas[1,1]
b2True <- trueBetas[2,1]
b3True <- trueBetas[3,1]

R2True <- 1 - sig2True/varY
```

#### Set Up - Input Simulation Parameters
```
library(tidyverse)

mu_x1 <- 1
mu_x2 <- 0.5
mu_x3 <- 1.5

numYears <- 3

# Sample size per study
#
n <- 50
#n <- 100
#n <- 150

if(case==1){
  numCountries <-5
  spreadMuY <- -2 #small
  
}else if(case==2){
  numCountries <-5
  spreadMuY <- -10 #large
  
}else if(case==3){
  numCountries <- 5
  spreadMuY <- -6 #medium
  
}else if(case==4){
  numCountries <-19
  spreadMuY <- -2 #small
  
}else if(case==5){
  numCountries <-19
  spreadMuY <- -10 #large
  
}else if(case==6){
  numCountries <-19
  spreadMuY <- -6 #medium
  
}else if(case==7){
  numCountries <-33
  spreadMuY <- -2 #small
  
}else if(case==8){
  numCountries <-33
  spreadMuY <- -10 #large
  
}else if(case==9){
  numCountries <-33
  spreadMuY <- -6 #small
}
# spreadMuY <- Will muY have a large variation (-10 to 10), or small variation (-2 to 2)?
```

#### Data Generation
```
library(tidyverse)
library(lmtest)
library(plm)
library(lme4)
library(sandwich)

#HOW THE DATA IS GENERATED
  #-Each (country) has 1 study ran per year for 5 years
  #-Heterogeneity in countries is captured by the different mean value for each country
  #-Heterogeneity in time is captured by the trend in y=y+0.5*i

###### COUNTRY 1 ######

### 2020 ###
mu_y <- spreadMuY
mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # The mean for y changes, x mean remains

z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
y <- z_data[, 1]
x <- z_data[, 2:4]

X <- data.frame(constant = 1,
                x = x, 
                year = 2020, 
                country = 1, 
                paper = 1,
                trend = -1)

### 2021-2024 ###
for (i in 1:(numYears-1)) {
  #mu_y <- mu_y+(timeHet*i) # The mean by _ per year
  #print(mu_y)
  mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x mean remains
  
  z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
  y2 <- z_data[, 1]
  x <- z_data[, 2:4]

  X2 <- data.frame(constant = 1,
                   x = x,
                   year = i + 2020,
                   country = 1,
                   paper = i + 1,
                   trend = -1+i/2)
  
  X <- bind_rows(X, X2)
  y <- c(y, y2)
}


  
###### COUNTRY 2 - 5 (or 9) ######

### 2020-2024 ###
paper_iterator <- numYears+1 # since Country 1 includes 1 paper per year
countryIndex=2

if(numCountries==5 & spreadMuY==-2){
  for (i in 2:numCountries) {
    for (j in 0:(numYears-1)) {
      mu_y <- i-3#+(timeHet*j) # The mean for y increases by 1 per country, also increases by _ per year
      #print(mu_y)
      
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x means remain
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = i,
                       paper = paper_iterator,
                       trend = -1+j/2)
      
      X <- bind_rows(X, X2)
      y <- c(y, y2)
      paper_iterator <- paper_iterator + 1
    }
  }
}
else if(numCountries==5 & spreadMuY==-10){
  for (i in seq(from=-5, to=10, by=5)) { # (-10), -5, ..., 10
    for (j in 0:(numYears-1)) {
      mu_y <- i#+(timeHet*j) # The mean for y increases by 5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==5 & spreadMuY==-6){
  for (i in seq(from=-3, to=6, by=3)) { # (-6), -3, ..., 6
    for (j in 0:(numYears-1)) {
      mu_y <- i#+(timeHet*j) # The mean for y increases by 5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}  

else if(numCountries==19 & spreadMuY==-2){
  temp=c(-1.76,-1.54,-1.32,-1.1,-0.88,-0.66,-0.44,-0.22,0,0.22,0.44,0.66,0.88,1.1,1.32,1.54,1.76,2) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==19 & spreadMuY==-10){
  temp=c(-8.8,-7.7,-6.6,-5.5,-4.4,-3.3,-2.2,-1.1,0,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,10) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==19 & spreadMuY==-6){
  temp=c(-5.28,-4.62,-3.96,-3.3,-2.64,-1.98,-1.32,-0.66,0,0.66,1.32,1.98,2.64,3.3,3.96,4.62,5.28,6) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==33 & spreadMuY==-2){
  temp=c(-1.875,-1.75,-1.625,-1.5,-1.375,-1.25,-1.125,-1,-0.875,-0.75,-0.625,-0.5,-0.375,-0.25,-0.125,0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1.125,1.25,1.375,1.5,1.625,1.75,1.875,2) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==33 & spreadMuY==-10){
  temp=c(-9.375,-8.75,-8.125,-7.5,-6.875,-6.25,-5.625,-5,-4.375,-3.75,-3.125,-2.5,-1.875,-1.25,-0.625,0,0.625,1.25,1.875,2.5,3.125,3.75,4.375,5,5.625,6.25,6.875,7.5,8.125,8.75,9.375,10) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}

else if(numCountries==33 & spreadMuY==-6){
  temp=c(-5.625,-5.25,-4.875,-4.5,-4.125,-3.75,-3.375,-3,-2.625,-2.25,-1.875,-1.5,-1.125,-0.75,-0.375,0,0.375,0.75,1.125,1.5,1.875,2.25,2.625,3,3.375,3.75,4.125,4.5,4.875,5.25,5.625,6) # for calculating spreadMuY
  for (i in 1:(length(temp))){
    for (j in 0:(numYears-1)) {
      mu_y <- temp[i]#+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
      #print(mu_y)
      mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2:4]
      y2 <- z_data[, 1]
      
      X2 <- data.frame(constant = 1,
                       x = x,
                       year = j + 2020,
                       country = countryIndex,
                       paper = paper_iterator,
                       trend = -1+j/2)
      paper_iterator <- paper_iterator + 1
      
      X <- rbind(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}
```

#### Run Models
```
# To ensure NAs don't keep the regressions from running
y <- na.omit(y)
X <- na.omit(X)

### Random Effects (Study Level) ###
RE <- plm(y ~ x.1 + x.2 + x.3, data=X, index=c("paper"), model="random")  #random model
#summary(RE)

REmse[k] <- mean(RE$residuals^2)                # MSE (df adjusted)
REmae[k] <- mean(abs(RE$residuals))             # MAE (df adjusted)
REmpe[k] <- mean(RE$residuals/RE$model$y)       # MPE (df adjusted)
REmape[k] <- abs(REmpe[k])

REb0[k] <- coef(RE)[1]
REb1[k] <- coef(RE)[2]
REb2[k] <- coef(RE)[3]
REb3[k] <- coef(RE)[4]

REmse_x1[k] <- mean((REb1[k]-b1True)^2)         # MSE (df adjusted)
REmae_x1[k] <- mean(abs(REb1[k]-b1True))        # MAE (df adjusted)
REmpe_x1[k] <- mean((b1True-REb1[k])/b1True)   # MPE (df adjusted)
REmape_x1[k] <- abs(REmpe_x1[k])
REmse_x2[k] <- mean((REb2[k]-b2True)^2)         # MSE (df adjusted)
REmae_x2[k] <- mean(abs(REb2[k]-b2True))        # MAE (df adjusted)
REmpe_x2[k] <- mean((b2True-REb2[k])/b2True)   # MPE (df adjusted)
REmape_x2[k] <- abs(REmpe_x2[k])
REmse_x3[k] <- mean((REb3[k]-b3True)^3)         # MSE (df adjusted)
REmae_x3[k] <- mean(abs(REb3[k]-b3True))        # MAE (df adjusted)
REmpe_x3[k] <- mean((b3True-REb3[k])/b3True)   # MPE (df adjusted)
REmape_x3[k] <- abs(REmpe_x3[k])

REsig2[k] <- sigma(RE)^2

REb0SE[k] <- sqrt(diag(vcov(RE)))[1]
REb1SE[k] <- sqrt(diag(vcov(RE)))[2]
REb2SE[k] <- sqrt(diag(vcov(RE)))[3]
REb3SE[k] <- sqrt(diag(vcov(RE)))[4]

CIlo_REb0[k] <- REb0[k] - (1.96*REb0SE[k])
CIhi_REb0[k] <- REb0[k] + (1.96*REb0SE[k])
CIlo_REb1[k] <- REb1[k] - (1.96*REb1SE[k])
CIhi_REb1[k] <- REb1[k] + (1.96*REb1SE[k])
CIlo_REb2[k] <- REb2[k] - (1.96*REb2SE[k])
CIhi_REb2[k] <- REb2[k] + (1.96*REb2SE[k])
CIlo_REb3[k] <- REb3[k] - (1.96*REb3SE[k])
CIhi_REb3[k] <- REb3[k] + (1.96*REb3SE[k])
    
### Study-level Fixed Effects ###
FE <- lm(y ~ x.1 + x.2 + x.3 + as.factor(paper), data=X)  #fixed model
FEresult <- summary(FE)

FEmse[k] <- mean(FE$residuals^2)                # MSE (df adjusted)
FEmae[k] <- mean(abs(FE$residuals))             # MAE (df adjusted)
FEmpe[k] <- mean(FE$residuals/FE$model$y)       # MPE (df adjusted)
FEmape[k] <- abs(FEmpe[k])

FEb0[k] <- coef(FE)[1]
FEb1[k] <- coef(FE)[2]
FEb2[k] <- coef(FE)[3]
FEb3[k] <- coef(FE)[4]

FEmse_x1[k] <- mean((FEb1[k]-b1True)^2)         # MSE (df adjusted)
FEmae_x1[k] <- mean(abs(FEb1[k]-b1True))        # MAE (df adjusted)
FEmpe_x1[k] <- mean((b1True-FEb1[k])/b1True)   # MPE (df adjusted)
FEmape_x1[k] <- abs(FEmpe_x1[k])
FEmse_x2[k] <- mean((FEb2[k]-b2True)^2)         # MSE (df adjusted)
FEmae_x2[k] <- mean(abs(FEb2[k]-b2True))        # MAE (df adjusted)
FEmpe_x2[k] <- mean((b2True-FEb2[k])/b2True)   # MPE (df adjusted)
FEmape_x2[k] <- abs(FEmpe_x2[k])
FEmse_x3[k] <- mean((FEb3[k]-b3True)^3)         # MSE (df adjusted)
FEmae_x3[k] <- mean(abs(FEb3[k]-b3True))        # MAE (df adjusted)
FEmpe_x3[k] <- mean((b3True-FEb3[k])/b3True)   # MPE (df adjusted)
FEmape_x3[k] <- abs(FEmpe_x3[k])

FEsig2[k] <- sigma(FE)^2

FEb0SE[k] <- sqrt(diag(vcov(FE)))[1]
FEb1SE[k] <- sqrt(diag(vcov(FE)))[2]
FEb2SE[k] <- sqrt(diag(vcov(FE)))[3]
FEb3SE[k] <- sqrt(diag(vcov(FE)))[4]

CIlo_FEb0[k] <- FEb0[k] - (1.96*FEb0SE[k])
CIhi_FEb0[k] <- FEb0[k] + (1.96*FEb0SE[k])
CIlo_FEb1[k] <- FEb1[k] - (1.96*FEb1SE[k])
CIhi_FEb1[k] <- FEb1[k] + (1.96*FEb1SE[k])
CIlo_FEb2[k] <- FEb2[k] - (1.96*FEb2SE[k])
CIhi_FEb2[k] <- FEb2[k] + (1.96*FEb2SE[k])
CIlo_FEb3[k] <- FEb3[k] - (1.96*FEb3SE[k])
CIhi_FEb3[k] <- FEb3[k] + (1.96*FEb3SE[k])
  
### Two Way Fixed Effects ###
FE2 <- lm(y ~ x.1 + x.2 + x.3 + as.factor(country) + as.factor(year), data = X)
FE2result <- summary(FE2)

FE2mse[k] <- mean(FE2$residuals^2)                # MSE (df adjusted)
FE2mae[k] <- mean(abs(FE2$residuals))             # MAE (df adjusted)
FE2mpe[k] <- mean(FE2$residuals/FE2$model$y)       # MPE (df adjusted)
FE2mape[k] <- abs(FE2mpe[k])

FE2b0[k] <- coef(FE2)[1]
FE2b1[k] <- coef(FE2)[2]
FE2b2[k] <- coef(FE2)[3]
FE2b3[k] <- coef(FE2)[4]

FE2mse_x1[k] <- mean((FE2b1[k]-b1True)^2)         # MSE (df adjusted)
FE2mae_x1[k] <- mean(abs(FE2b1[k]-b1True))        # MAE (df adjusted)
FE2mpe_x1[k] <- mean((b1True-FE2b1[k])/b1True)   # MPE (df adjusted)
FE2mape_x1[k] <- abs(FE2mpe_x1[k])
FE2mse_x2[k] <- mean((FE2b2[k]-b2True)^2)         # MSE (df adjusted)
FE2mae_x2[k] <- mean(abs(FE2b2[k]-b2True))        # MAE (df adjusted)
FE2mpe_x2[k] <- mean((b2True-FE2b2[k])/b2True)   # MPE (df adjusted)
FE2mape_x2[k] <- abs(FE2mpe_x2[k])
FE2mse_x3[k] <- mean((FE2b3[k]-b3True)^3)         # MSE (df adjusted)
FE2mae_x3[k] <- mean(abs(FE2b3[k]-b3True))        # MAE (df adjusted)
FE2mpe_x3[k] <- mean((b3True-FE2b3[k])/b3True)   # MPE (df adjusted)
FE2mape_x3[k] <- abs(FE2mpe_x3[k])

FE2sig2[k] <- sigma(FE2)^2

FE2b0SE[k] <- sqrt(diag(vcov(FE2)))[1]
FE2b1SE[k] <- sqrt(diag(vcov(FE2)))[2]
FE2b2SE[k] <- sqrt(diag(vcov(FE2)))[3]
FE2b3SE[k] <- sqrt(diag(vcov(FE2)))[4]

CIlo_FE2b0[k] <- FE2b0[k] - (1.96*FE2b0SE[k])
CIhi_FE2b0[k] <- FE2b0[k] + (1.96*FE2b0SE[k])
CIlo_FE2b1[k] <- FE2b1[k] - (1.96*FE2b1SE[k])
CIhi_FE2b1[k] <- FE2b1[k] + (1.96*FE2b1SE[k])
CIlo_FE2b2[k] <- FE2b2[k] - (1.96*FE2b2SE[k])
CIhi_FE2b2[k] <- FE2b2[k] + (1.96*FE2b2SE[k])
CIlo_FE2b3[k] <- FE2b3[k] - (1.96*FE2b3SE[k])
CIhi_FE2b3[k] <- FE2b3[k] + (1.96*FE2b3SE[k])
  
### [LOCATION-level] Random Effects ###
REl <- plm(y ~ x.1 + x.2 + x.3, data=X, index=c("country"), model="random")  #random model
#summary(RE)

REl_mse[k] <- mean(REl$residuals^2)                # MSE (df adjusted)
REl_mae[k] <- mean(abs(REl$residuals))             # MAE (df adjusted)
REl_mpe[k] <- mean(REl$residuals/REl$model$y)       # MPE (df adjusted)
REl_mape[k] <- abs(REl_mpe[k])

REl_b0[k] <- coef(REl)[1]
REl_b1[k] <- coef(REl)[2]
REl_b2[k] <- coef(REl)[3]
REl_b3[k] <- coef(REl)[4]

REl_mse_x1[k] <- mean((REl_b1[k]-b1True)^2)         # MSE (df adjusted)
REl_mae_x1[k] <- mean(abs(REl_b1[k]-b1True))        # MAE (df adjusted)
REl_mpe_x1[k] <- mean((b1True-REl_b1[k])/b1True)   # MPE (df adjusted)
REl_mape_x1[k] <- abs(REl_mpe_x1[k])
REl_mse_x2[k] <- mean((REl_b2[k]-b2True)^2)         # MSE (df adjusted)
REl_mae_x2[k] <- mean(abs(REl_b2[k]-b2True))        # MAE (df adjusted)
REl_mpe_x2[k] <- mean((b2True-REl_b2[k])/b2True)   # MPE (df adjusted)
REl_mape_x2[k] <- abs(REl_mpe_x2[k])
REl_mse_x3[k] <- mean((REl_b3[k]-b3True)^3)         # MSE (df adjusted)
REl_mae_x3[k] <- mean(abs(REl_b3[k]-b3True))        # MAE (df adjusted)
REl_mpe_x3[k] <- mean((b3True-REl_b3[k])/b3True)   # MPE (df adjusted)
REl_mape_x3[k] <- abs(REl_mpe_x3[k])

REl_sig2[k] <- sigma(REl)^2

REl_b0SE[k] <- sqrt(diag(vcov(REl)))[1]
REl_b1SE[k] <- sqrt(diag(vcov(REl)))[2]
REl_b2SE[k] <- sqrt(diag(vcov(REl)))[3]
REl_b3SE[k] <- sqrt(diag(vcov(REl)))[4]

CIlo_REl_b0[k] <- REl_b0[k] - (1.96*REl_b0SE[k])
CIhi_REl_b0[k] <- REl_b0[k] + (1.96*REl_b0SE[k])
CIlo_REl_b1[k] <- REl_b1[k] - (1.96*REl_b1SE[k])
CIhi_REl_b1[k] <- REl_b1[k] + (1.96*REl_b1SE[k])
CIlo_REl_b2[k] <- REl_b2[k] - (1.96*REl_b2SE[k])
CIhi_REl_b2[k] <- REl_b2[k] + (1.96*REl_b2SE[k])
CIlo_REl_b3[k] <- REl_b3[k] - (1.96*REl_b3SE[k])
CIhi_REl_b3[k] <- REl_b3[k] + (1.96*REl_b3SE[k])
  
### [LOCATION-level] Fixed Effects ###
FEl <- lm(y ~ x.1 + x.2 + x.3 + as.factor(country), data=X)  #fixed model
FElresult <- summary(FEl)

FEl_mse[k] <- mean(FEl$residuals^2)                # MSE (df adjusted)
FEl_mae[k] <- mean(abs(FEl$residuals))             # MAE (df adjusted)
FEl_mpe[k] <- mean(FEl$residuals/FEl$model$y)       # MPE (df adjusted)
FEl_mape[k] <- abs(FEl_mpe[k])

FEl_b0[k] <- coef(FEl)[1]
FEl_b1[k] <- coef(FEl)[2]
FEl_b2[k] <- coef(FEl)[3]
FEl_b3[k] <- coef(FEl)[4]

FEl_mse_x1[k] <- mean((FEl_b1[k]-b1True)^2)         # MSE (df adjusted)
FEl_mae_x1[k] <- mean(abs(FEl_b1[k]-b1True))        # MAE (df adjusted)
FEl_mpe_x1[k] <- mean((b1True-FEl_b1[k])/b1True)   # MPE (df adjusted)
FEl_mape_x1[k] <- abs(FEl_mpe_x1[k])
FEl_mse_x2[k] <- mean((FEl_b2[k]-b2True)^2)         # MSE (df adjusted)
FEl_mae_x2[k] <- mean(abs(FEl_b2[k]-b2True))        # MAE (df adjusted)
FEl_mpe_x2[k] <- mean((b2True-FEl_b2[k])/b2True)   # MPE (df adjusted)
FEl_mape_x2[k] <- abs(FEl_mpe_x2[k])
FEl_mse_x3[k] <- mean((FEl_b3[k]-b3True)^3)         # MSE (df adjusted)
FEl_mae_x3[k] <- mean(abs(FEl_b3[k]-b3True))        # MAE (df adjusted)
FEl_mpe_x3[k] <- mean((b3True-FEl_b3[k])/b3True)   # MPE (df adjusted)
FEl_mape_x3[k] <- abs(FEl_mpe_x3[k])

FEl_sig2[k] <- sigma(FEl)^2

FEl_b0SE[k] <- sqrt(diag(vcov(FEl)))[1]
FEl_b1SE[k] <- sqrt(diag(vcov(FEl)))[2]
FEl_b2SE[k] <- sqrt(diag(vcov(FEl)))[3]
FEl_b3SE[k] <- sqrt(diag(vcov(FEl)))[4]

CIlo_FEl_b0[k] <- FEl_b0[k] - (1.96*FEl_b0SE[k])
CIhi_FEl_b0[k] <- FEl_b0[k] + (1.96*FEl_b0SE[k])
CIlo_FEl_b1[k] <- FEl_b1[k] - (1.96*FEl_b1SE[k])
CIhi_FEl_b1[k] <- FEl_b1[k] + (1.96*FEl_b1SE[k])
CIlo_FEl_b2[k] <- FEl_b2[k] - (1.96*FEl_b2SE[k])
CIhi_FEl_b2[k] <- FEl_b2[k] + (1.96*FEl_b2SE[k])
CIlo_FEl_b3[k] <- FEl_b3[k] - (1.96*FEl_b3SE[k])
CIhi_FEl_b3[k] <- FEl_b3[k] + (1.96*FEl_b3SE[k])  

### [Study-level] Fixed Effects (with Trend) ###
FEsT <- lm(y ~ x.1 + x.2 + x.3 + trend + as.factor(paper), data=X)  #fixed model
FEsTresult <- summary(FEsT)

FEsT_mse[k] <- mean(FEsT$residuals^2)                # MSE (df adjusted)
FEsT_mae[k] <- mean(abs(FEsT$residuals))             # MAE (df adjusted)
FEsT_mpe[k] <- mean(FEsT$residuals/FEsT$model$y)       # MPE (df adjusted)
FEsT_mape[k] <- abs(FEsT_mpe[k])

FEsT_b0[k] <- coef(FEsT)[1]
FEsT_b1[k] <- coef(FEsT)[2]
FEsT_b2[k] <- coef(FEsT)[3]
FEsT_b3[k] <- coef(FEsT)[4]
FEsT_bTrend[k] <- coef(FEsT)[5]

FEsT_mse_x1[k] <- mean((FEsT_b1[k]-b1True)^2)         # MSE (df adjusted)
FEsT_mae_x1[k] <- mean(abs(FEsT_b1[k]-b1True))        # MAE (df adjusted)
FEsT_mpe_x1[k] <- mean((b1True-FEsT_b1[k])/b1True)   # MPE (df adjusted)
FEsT_mape_x1[k] <- abs(FEsT_mpe_x1[k])
FEsT_mse_trend[k] <- mean((FEsT_bTrend[k]-0)^2)         # MSE (df adjusted)
FEsT_mae_trend[k] <- mean(abs(FEsT_bTrend[k]-0))        # MAE (df adjusted)
FEsT_mpe_trend[k] <- 0#mean((timeHet-FEsT_bTrend[k])/timeHet)   # MPE (df adjusted)
FEsT_mape_trend[k] <- 0#abs(FEsT_mpe_trend[k])
FEsT_mse_x2[k] <- mean((FEsT_b2[k]-b2True)^2)         # MSE (df adjusted)
FEsT_mae_x2[k] <- mean(abs(FEsT_b2[k]-b2True))        # MAE (df adjusted)
FEsT_mpe_x2[k] <- mean((b2True-FEsT_b2[k])/b2True)   # MPE (df adjusted)
FEsT_mape_x2[k] <- abs(FEsT_mpe_x2[k])
FEsT_mse_x3[k] <- mean((FEsT_b3[k]-b3True)^3)         # MSE (df adjusted)
FEsT_mae_x3[k] <- mean(abs(FEsT_b3[k]-b3True))        # MAE (df adjusted)
FEsT_mpe_x3[k] <- mean((b3True-FEsT_b3[k])/b3True)   # MPE (df adjusted)
FEsT_mape_x3[k] <- abs(FEsT_mpe_x3[k])

FEsT_sig2[k] <- sigma(FEsT)^2

FEsT_b0SE[k] <- sqrt(diag(vcov(FEsT)))[1]
FEsT_b1SE[k] <- sqrt(diag(vcov(FEsT)))[2]
FEsT_b2SE[k] <- sqrt(diag(vcov(FEsT)))[3]
FEsT_b3SE[k] <- sqrt(diag(vcov(FEsT)))[4]
FEsT_bTrendSE[k] <- sqrt(diag(vcov(FEsT)))[5]

CIlo_FEsT_b0[k] <- FEsT_b0[k] - (1.96*FEsT_b0SE[k])
CIhi_FEsT_b0[k] <- FEsT_b0[k] + (1.96*FEsT_b0SE[k])
CIlo_FEsT_b1[k] <- FEsT_b1[k] - (1.96*FEsT_b1SE[k])
CIhi_FEsT_b1[k] <- FEsT_b1[k] + (1.96*FEsT_b1SE[k])
CIlo_FEsT_b2[k] <- FEsT_b2[k] - (1.96*FEsT_b2SE[k])
CIhi_FEsT_b2[k] <- FEsT_b2[k] + (1.96*FEsT_b2SE[k])
CIlo_FEsT_b3[k] <- FEsT_b3[k] - (1.96*FEsT_b3SE[k])
CIhi_FEsT_b3[k] <- FEsT_b3[k] + (1.96*FEsT_b3SE[k])
CIlo_FEsT_bTrend[k] <- FEsT_bTrend[k] - (1.96*FEsT_bTrendSE[k])
CIhi_FEsT_bTrend[k] <- FEsT_bTrend[k] + (1.96*FEsT_bTrendSE[k])
      
### [Location-level] Fixed Effects (with Trend) ###
FElT <- lm(y ~ x.1 + x.2 + x.3 + trend + as.factor(country), data=X)  #fixed model
FElTresult <- summary(FElT)

FElT_mse[k] <- mean(FElT$residuals^2)                # MSE (df adjusted)
FElT_mae[k] <- mean(abs(FElT$residuals))             # MAE (df adjusted)
FElT_mpe[k] <- mean(FElT$residuals/FElT$model$y)       # MPE (df adjusted)
FElT_mape[k] <- abs(FElT_mpe[k])

FElT_b0[k] <- coef(FElT)[1]
FElT_b1[k] <- coef(FElT)[2]
FElT_b2[k] <- coef(FElT)[3]
FElT_b3[k] <- coef(FElT)[4]
FElT_bTrend[k] <- coef(FElT)[5]

FElT_mse_x1[k] <- mean((FElT_b1[k]-b1True)^2)         # MSE (df adjusted)
FElT_mae_x1[k] <- mean(abs(FElT_b1[k]-b1True))        # MAE (df adjusted)
FElT_mpe_x1[k] <- mean((b1True-FElT_b1[k])/b1True)   # MPE (df adjusted)
FElT_mape_x1[k] <- abs(FElT_mpe_x1[k])
FElT_mse_trend[k] <- mean((FElT_bTrend[k]-0)^2)         # MSE (df adjusted)
FElT_mae_trend[k] <- mean(abs(FElT_bTrend[k]-0))        # MAE (df adjusted)
FElT_mpe_trend[k] <- 0#mean((timeHet-FElT_bTrend[k])/timeHet)   # MPE (df adjusted)
FElT_mape_trend[k] <- 0#abs(FElT_mpe_trend[k])
FElT_mse_x2[k] <- mean((FElT_b2[k]-b2True)^2)         # MSE (df adjusted)
FElT_mae_x2[k] <- mean(abs(FElT_b2[k]-b2True))        # MAE (df adjusted)
FElT_mpe_x2[k] <- mean((b2True-FElT_b2[k])/b2True)   # MPE (df adjusted)
FElT_mape_x2[k] <- abs(FElT_mpe_x2[k])
FElT_mse_x3[k] <- mean((FElT_b3[k]-b3True)^3)         # MSE (df adjusted)
FElT_mae_x3[k] <- mean(abs(FElT_b3[k]-b3True))        # MAE (df adjusted)
FElT_mpe_x3[k] <- mean((b3True-FElT_b3[k])/b3True)   # MPE (df adjusted)
FElT_mape_x3[k] <- abs(FElT_mpe_x3[k])

FElT_sig2[k] <- sigma(FElT)^2

FElT_b0SE[k] <- sqrt(diag(vcov(FElT)))[1]
FElT_b1SE[k] <- sqrt(diag(vcov(FElT)))[2]
FElT_b2SE[k] <- sqrt(diag(vcov(FElT)))[3]
FElT_b3SE[k] <- sqrt(diag(vcov(FElT)))[4]
FElT_bTrendSE[k] <- sqrt(diag(vcov(FElT)))[5]

CIlo_FElT_b0[k] <- FElT_b0[k] - (1.96*FElT_b0SE[k])
CIhi_FElT_b0[k] <- FElT_b0[k] + (1.96*FElT_b0SE[k])
CIlo_FElT_b1[k] <- FElT_b1[k] - (1.96*FElT_b1SE[k])
CIhi_FElT_b1[k] <- FElT_b1[k] + (1.96*FElT_b1SE[k])
CIlo_FElT_b2[k] <- FElT_b2[k] - (1.96*FElT_b2SE[k])
CIhi_FElT_b2[k] <- FElT_b2[k] + (1.96*FElT_b2SE[k])
CIlo_FElT_b3[k] <- FElT_b3[k] - (1.96*FElT_b3SE[k])
CIhi_FElT_b3[k] <- FElT_b3[k] + (1.96*FElT_b3SE[k])
CIlo_FElT_bTrend[k] <- FElT_bTrend[k] - (1.96*FElT_bTrendSE[k])
CIhi_FElT_bTrend[k] <- FElT_bTrend[k] + (1.96*FElT_bTrendSE[k])
  
### Mixed Effects (RE @ Study-level) ###
ME <- lmer(y ~ x.1 + x.2 + x.3 + (1 | paper), data=X)  #fixed model
MEresult <- summary(ME)

MEmse[k] <- mean(residuals(ME)^2)                # MSE (df adjusted)
MEmae[k] <- mean(abs(residuals(ME)))             # MAE (df adjusted)
MEmpe[k] <- mean(residuals(ME)/y)            # MPE (df adjusted)
MEmape[k] <- abs(MEmpe[k])

MEb0[k] <- summary(ME)$coefficients[1]
MEb1[k] <- summary(ME)$coefficients[2]
MEb2[k] <- summary(ME)$coefficients[3]
MEb3[k] <- summary(ME)$coefficients[4]

ME_mse_x1[k] <- mean((MEb1[k]-b1True)^2)         # MSE (df adjusted)
ME_mae_x1[k] <- mean(abs(MEb1[k]-b1True))        # MAE (df adjusted)
ME_mpe_x1[k] <- mean((b1True-MEb1[k])/b1True)   # MPE (df adjusted)
ME_mape_x1[k] <- abs(ME_mpe_x1[k])
ME_mse_x2[k] <- mean((MEb2[k]-b2True)^2)         # MSE (df adjusted)
ME_mae_x2[k] <- mean(abs(MEb2[k]-b2True))        # MAE (df adjusted)
ME_mpe_x2[k] <- mean((b2True-MEb2[k])/b2True)   # MPE (df adjusted)
ME_mape_x2[k] <- abs(ME_mpe_x2[k])
ME_mse_x3[k] <- mean((MEb3[k]-b3True)^3)         # MSE (df adjusted)
ME_mae_x3[k] <- mean(abs(MEb3[k]-b3True))        # MAE (df adjusted)
ME_mpe_x3[k] <- mean((b3True-MEb3[k])/b3True)   # MPE (df adjusted)
ME_mape_x3[k] <- abs(ME_mpe_x3[k])

MEsig2[k] <- sigma(ME)^2

MEb0SE[k] <- summary(ME)$coefficients[5]
MEb1SE[k] <- summary(ME)$coefficients[6]
MEb2SE[k] <- summary(ME)$coefficients[7]
MEb3SE[k] <- summary(ME)$coefficients[8]

CIlo_MEb0[k] <- MEb0[k] - (1.96*MEb0SE[k])
CIhi_MEb0[k] <- MEb0[k] + (1.96*MEb0SE[k])
CIlo_MEb1[k] <- MEb1[k] - (1.96*MEb1SE[k])
CIhi_MEb1[k] <- MEb1[k] + (1.96*MEb1SE[k])
CIlo_MEb2[k] <- MEb2[k] - (1.96*MEb2SE[k])
CIhi_MEb2[k] <- MEb2[k] + (1.96*MEb2SE[k])
CIlo_MEb3[k] <- MEb3[k] - (1.96*MEb3SE[k])
CIhi_MEb3[k] <- MEb3[k] + (1.96*MEb3SE[k])
    
### Mixed Effects (RE @ Country-level) ###
MEl <- lmer(y ~ x.1 + x.2 + x.3 + (1 | country), data=X)  #fixed model
MElresult <- summary(MEl)

MEl_mse[k] <- mean(residuals(MEl)^2)                        # MSE (df adjusted)
MEl_mae[k] <- mean(abs(residuals(MEl)))             # MAE (df adjusted)
MEl_mpe[k] <- mean(residuals(MEl)/y)            # MPE (df adjusted)
MEl_mape[k] <- abs(MEl_mpe[k])

MEl_b0[k] <- summary(MEl)$coefficients[1]
MEl_b1[k] <- summary(MEl)$coefficients[2]
MEl_b2[k] <- summary(MEl)$coefficients[3]
MEl_b3[k] <- summary(MEl)$coefficients[4]

MEl_mse_x1[k] <- mean((MEl_b1[k]-b1True)^2)         # MSE (df adjusted)
MEl_mae_x1[k] <- mean(abs(MEl_b1[k]-b1True))        # MAE (df adjusted)
MEl_mpe_x1[k] <- mean((b1True-MEl_b1[k])/b1True)   # MPE (df adjusted)
MEl_mape_x1[k] <- abs(MEl_mpe_x1[k])
MEl_mse_x2[k] <- mean((MEl_b2[k]-b2True)^2)         # MSE (df adjusted)
MEl_mae_x2[k] <- mean(abs(MEl_b2[k]-b2True))        # MAE (df adjusted)
MEl_mpe_x2[k] <- mean((b2True-MEl_b2[k])/b2True)   # MPE (df adjusted)
MEl_mape_x2[k] <- abs(MEl_mpe_x2[k])
MEl_mse_x3[k] <- mean((MEl_b3[k]-b3True)^3)         # MSE (df adjusted)
MEl_mae_x3[k] <- mean(abs(MEl_b3[k]-b3True))        # MAE (df adjusted)
MEl_mpe_x3[k] <- mean((b3True-MEl_b3[k])/b3True)   # MPE (df adjusted)
MEl_mape_x3[k] <- abs(MEl_mpe_x3[k])

MEl_sig2[k] <- sigma(MEl)^2

MEl_b0SE[k] <- summary(MEl)$coefficients[5]
MEl_b1SE[k] <- summary(MEl)$coefficients[6]
MEl_b2SE[k] <- summary(MEl)$coefficients[7]
MEl_b3SE[k] <- summary(MEl)$coefficients[8]

CIlo_MEl_b0[k] <- MEl_b0[k] - (1.96*MEl_b0SE[k])
CIhi_MEl_b0[k] <- MEl_b0[k] + (1.96*MEl_b0SE[k])
CIlo_MEl_b1[k] <- MEl_b1[k] - (1.96*MEl_b1SE[k])
CIhi_MEl_b1[k] <- MEl_b1[k] + (1.96*MEl_b1SE[k])
CIlo_MEl_b2[k] <- MEl_b2[k] - (1.96*MEl_b2SE[k])
CIhi_MEl_b2[k] <- MEl_b2[k] + (1.96*MEl_b2SE[k])
CIlo_MEl_b3[k] <- MEl_b3[k] - (1.96*MEl_b3SE[k])
CIhi_MEl_b3[k] <- MEl_b3[k] + (1.96*MEl_b3SE[k])
```
