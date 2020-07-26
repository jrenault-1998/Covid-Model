# clear workspace
rm(list=ls())
# this is the R package to solve the system of ODEs 
require(deSolve)
# need to invert matrices
require(matlib)
library(data.table)
library(curl)

# Set this working directory to where you want to save the figures.
# You need to write your own path here
###setwd('~/Desktop/MUN COVID/Second Wave/Figures/')

# The code below is commented out because it's just an illustrative graph
# about something else
# x = seq(1,5,1)
# x2 = seq(1,5,.1)
# y = 10*exp(-(x-1))
# par(mar=c(5,3,2,2)+0.1)
# png(file="well-being.png",width=500,height=450)
# plot(x,y,pch=19, xlab = "Alert level", ylab = "", bty="n",yaxt="n", cex.lab=2, cex = 1.5, cex.axis=1.5)
# title(ylab="Well-being", line=1, cex.lab=2)
# axis(side=2,labels = FALSE, tck=0)
# lines(x2,10*exp(-(x2-1)),lwd=1.5)
# lines(c(1,5),c(10,10*exp(-4)),lty=2,lwd=2)
# dev.off()

#data <- data.table::fread('https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/COVID19_Canada.csv', fill = TRUE)[Province == 'ON']

# Parameters which should be a direct copy of Miller et al.
# See Miller et al. for definitions.
bA <- 0.5
deltaE <- 1/4
deltaP <- 1/2.4
deltaC <- 1/3.2
deltaA <- 1/7
# Rate of quarantine due to contact tracing
lambda = 1/5
# R0 after re-escalation
R0c <- 0.5
bC <- 0.1
N<-1

# Added importations since initial conditions of one exposed gives
# very slow initial growth of the outbreak.
tau <- 50/519705

# vector of age-structured asymptomatic rates from Miller et al.
ri <- c(0.056, 0.056, 0.129, 0.281, 0.429, 0.490, 0.740, 0.740, 0.740)

# vector of age-structured hospitalization rates for symptomatic cases
h_r = c(0.001, 0.003, 0.012, 0.032, 0.049, 0.102, 0.166, 0.243, 0.273)

# vector of age-structured ICU rates given hospitalization
i_r = c(0.05, 0.05, 0.05, 0.05, 0.063, 0.122, 0.274, 0.432, 0.709)

# Number of NLers in each age category from 2016 census data
n_09 = 22365+26040
n_1019 = 26035+27255
n_2029 = 27700+28245
n_3039 = 29405+30740
n_4049 = 34505+38665
n_5059 = 42620+43080
n_6069 = 42025+37485
n_7079 = 26170+16950
n_80p = 11060+9360

# Vectorize and convert to fractions
age.frac = c(n_09, n_1019, n_2029, n_3039, n_4049, n_5059, n_6069, n_7079,n_80p)
tot.popn = sum(age.frac)
age.frac = age.frac/tot.popn

# Calculate average r, hospitalization given clinical infection, ICU given hospitalization
# Early on, I had been considering doing this model with age-structure, but I ended up
# deciding that the age-structure wasn't a critical piece, so it gets removed here.
r <- sum(ri*age.frac)
# Overwrite the calculation of r to be consistent with JC's work - specifically,
# we found that for the NL data, if it were assumed that asymptomatics were
# greater than 20% of infections the fit of this model to data was poor.
r <- 0.85
av.hosp <- sum(h_r*age.frac)
av.ICU <-  sum(i_r*age.frac)

# Given a user supplied R0, calculate the beta
find.beta = function(beta,R0){
  # Calculate R0 using the next-generation approach of van den Driesshe and
  # Watmough 2002.
F <- matrix(rep(0,16), nrow =4, ncol = 4)
V <- F
F[1,] <- c(0, beta, beta*bC, beta*bA)
V[1,1] <- deltaE + lambda
V[2,1] <- -r*deltaE
V[2,2] <- deltaP + lambda
V[3,2] <- -deltaP
V[3,3] <- deltaC + lambda
V[4,1] <- -(1-r)*deltaE
V[4,4] <- deltaA + lambda
inv.V = inv(V)
K = F%*%inv.V
R0.calc = max(abs(eigen(K)$values))
R0 - R0.calc
}

# Root finding to get the beta corresponding to the target R0
# Some of the code in this section is copy & paste from the Miller et al. code

betac <- uniroot(find.beta,c(0,10),R0=R0c)$root
print(betac)

# This is a function describing importations. This states that the importation
# rate is turned off on day 5. Without this forcing the epidemic begins too slowly.
tau.fun = function(t){
  if(t>=5){
    tauval <- 0
  } else{
    tauval <- tau}
  return(tauval)
}

# In NL, once stay-at-home orders were issued the outbreak was brought under
# control. Therefore, the transmission rate, beta, is assumed to have a breakpoint
# around t = SDstart which is a time at which stay-at-home orders are issued.
beta.fun = function(t){
  if(t>=SDstart){
    betaval <- betac
  } else{
  betaval <- beta}
  return(betaval)
}

# This is the system of ODEs for Miller et al., but modified so that individuals
# go into isolation due to contact tracing at rate lambda.
Miller.CT = function(t,y,parms){
  S <- y[1]
  E <- y[2]
  IP <- y[3]
  IC <- y[4]
  IA <- y[5]
  RS <- y[6]
  RA <- y[7]
  
  dS = -S*beta.fun(t)*(IP+bC*IC+bA*IA)/N
  dE = tau.fun(t)+S*beta.fun(t)*(IP+bC*IC+bA*IA)/N-deltaE*E - lambda*E
  dIP = r*deltaE*E-deltaP*IP - lambda*IP
  dIC = deltaP*IP-deltaC*IC - lambda*IC
  dIA = (1-r)*deltaE*E - deltaA*IA - lambda*IA
  dRS = deltaC*IC - r*lambda*E + lambda*IP + lambda*IC
  dRA = deltaA*IA - (1-r)*lambda*E + lambda*IA
  # cumulative clinical cases
  dcumIC = deltaP*IP
  
  return(list(c(dS,dE,dIP,dIC,dIA,dRS,dRA,dcumIC)))
}

# Assumed the epidemic begins with 3 exponsed.
E0 = 3/tot.popn
yini  = c(S = 1, E = E0, IP = 0, IC = 0, IA = 0, RS=0, RA=0, cumIC = 0)
# the times for the numerical integration
times <- seq(0, 365, by = .1)
# The assumed R0 prior to the higher alert level
R0 <-2.5
beta<-uniroot(find.beta,c(0,10),R0=R0)$root

# Day 7 for the implementation of the higher alert level
SDstart <- 7
out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
out7 <- data.frame(out)

# Day 10 for the implementation of the higher alert level
SDstart <- 10
out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
out10 <- data.frame(out)

# Making the figure for escalation on day 7 or day 10
png(file="reescalation.png",width=1000,height=450)
par(mfrow = c(1,2), mar=c(5,6,3,1))
# New daily cases is the derivative of cumulative cases
plot(tail(out7$time,-1), diff(out7$cumIC)*tot.popn, typ = "l", xlab = "days", ylab="daily cases", xlim = c(0,30), main = "Re-escalation on day 7", ylim = c(0,15), bty="n", lwd =3,cex.lab=2, cex.main = 2, cex.axis=2)
lines(c(7, 7), c(0,max(out7$IC*tot.popn)), lty = 2, lwd=2)
plot(tail(out10$time,-1), diff(out10$cumIC*tot.popn), typ = "l", xlab = "days", ylab="daily cases", xlim = c(0,30), main = "Re-escalation on day 10",ylim = c(0,15), bty="n", lwd=3,cex.lab=2, cex.main = 2, cex.axis=2)
lines(c(10, 10), c(0,max(out10$IC*tot.popn)), lty = 2, lwd=2)
dev.off()

# Loop across escalation days
SDstart.vec <- seq(0,20,1)
# Preallocate output vector
out.dat = data.frame("Estart"=NULL)
Estart = NULL
ICstart = NULL
Efinal = NULL
ICfinal = NULL
LDend = NULL

# The is a loop across different days that escalation could start
for(i in seq(1,length(SDstart.vec))){
  SDstart <- SDstart.vec[i]
  # performing the numerical integration
  out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
  out <- data.frame(out)
  # recording the number of exposed when re-escalation starts
  Estart[i] = out$E[which(out$time == SDstart)]
  # recording the number of clinically infected when re-escalation starts
  ICstart[i] = out$IC[which(out$time == SDstart)]
  # recording the number of exposed & clinically infected at the final time
  Efinal[i] = tail(out$cumE,1)
  ICfinal[i] = tail(out$cumIC,1)
  # the lockdown is assumed to end when the number of clinically infected
  # is < 1. The duration of the lockdown is the end time - the start time.
  LDend[i] = out$time[max(which(out$IC*tot.popn>1))]-SDstart
}

# Figure code
png(file="dontwait.png",width=1000,height=450)
par(mfrow = c(1,2), mar=c(5,6,3,1))
plot(SDstart.vec, tot.popn*ICfinal, ylab = "total clinical infections", xlab = "Days waited till re-escalation start", typ="l", ylim = c(0,3000),bty="n", lwd =3,cex.lab=2, cex.main = 2, cex.axis=2, main = "Delays are bad for public health")
plot(SDstart.vec, LDend, typ = "l", ylab = "Days at higher alert level", xlab = "Days waited till re-escalation start", bty="n", lwd =3,cex.lab=2, cex.main = 2, cex.axis=2, main = "Delays are bad for well-being")
dev.off()
L = length(SDstart)
mod = lm(LDend[5:L]~SDstart.vec[5:L])
# Slope is 0.41

# Below the code investigates different R0 values
R0vec<-c(1.3, 1.2, 1.1, 1.05)
LDend = matrix(rep(0,length(SDstart.vec)*length(R0vec)),nrow = length(SDstart.vec),ncol = length(R0vec))
for(j in seq(1,length(R0vec))){
  R0 = R0vec[j]
  beta<-uniroot(find.beta,c(0,10),R0=R0)$root
for(i in seq(1,length(SDstart.vec))){
  SDstart <- SDstart.vec[i]
  # performing the numerical integration
  out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
  out <- data.frame(out)
  LDend[i,j] = out$time[max(which(out$IC*tot.popn>1))]-SDstart
}
}

png(file="lowR0.png",width=500,height=450)
plot(SDstart.vec, LDend[,1], typ = "l", ylab = "lockdown duration (days)", xlab = "Days waited till lockdown start", ylim = c(28,38))
lines(SDstart.vec, LDend[,2])
lines(SDstart.vec, LDend[,3])
dev.off()

# The code below should be ignored. This is a start on writing some code for
# contact tracing with limited resources such that when too many individuals
# are infected contact tracing efficacy declines. The code below has errors, however,
# I think.
# Lmax = function(t){
#   if(t > t.on & t < t.off ){
#     Lmax = Lmax.val
#   } else {
#     Lmax = 0
#   }
#   return(Lmax)
# }

# Miller.CT.cap = function(t,y,parms){
#   S <- y[1]
#   E <- y[2]
#   IP <- y[3]
#   IC <- y[4]
#   IA <- y[5]
#   RS <- y[6]
#   RA <- y[7]
#   cumIC <-y[8]
#   dS = -S*beta.fun(t)*(IP+bC*IC+bA*IA)/N
#   if(tau.fun(t)+S*beta.fun(t)*(IP+bC*IC+bA*IA)/N <= Lmax(t)){
#     # Take everyone who would have been in E and put them in Q
#     dE = -deltaE*E
#   } else {
#     dE = (tau.fun(t)+S*beta.fun(t)*(IP+bC*IC+bA*IA)/N)*((tau.fun(t)+S*beta.fun(t)*(IP+bC*IC+bA*IA)/N)-Lmax(t))/(tau.fun(t)+S*beta.fun(t)*(IP+bC*IC+bA*IA)/N)-deltaE*E
#   }
#   dIP = r*deltaE*E-deltaP*IP
#   dIC = deltaP*IP-deltaC*IC
#   dIA = (1-r)*deltaE*E - deltaA*IA
#   dRS = deltaC*IC
#   dRA = deltaA*IA
#   # cumulative clinical cases
#   dcumIC = deltaP*IP
#   
#   return(list(c(dS,dE,dIP,dIC,dIA,dRS,dRA,dcumIC)))
# }

# 
# R0=2
# beta<-uniroot(find.beta,c(0,10),R0=R0)$root
# SDStart<-10
# 
# t.on <- 1
# t.off <- 15
# Lmax.val = (30-t.off+t.on)/tot.popn
# times = seq(0,100)
# yini  = c(S = 1, E = E0, IP = 0, IC = 0, IA = 0, RS=0, RA=0, cumIC=0)
# out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
# out <- data.frame(out)
# plot(out$time, out$IC*tot.popn, typ = "l", ylim = c(0,20), col = "red")
# 
# t.on <- 1
# t.off <- 7
# Lmax.val = (30-t.off+t.on)/tot.popn
# out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
# out <- data.frame(out)
# lines(out$time, out$IC*tot.popn, typ = "l", col="red", lty=2)
# lines(c(0,100), c(1,1))
# 
# t.on <- 8
# t.off <- 22
# Lmax.val = (30-t.off+t.on)/tot.popn
# out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
# out <- data.frame(out)
# lines(out$time, out$IC*tot.popn, typ = "l", col="orange")
# lines(c(0,100), c(1,1))
# 
# t.on <- 15
# t.off <- 29
# Lmax.val = (30-t.off+t.on)/tot.popn
# out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
# out <- data.frame(out)
# lines(out$time, out$IC*tot.popn, typ = "l", col="blue")
# 
# 

