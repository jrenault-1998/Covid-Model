# clear workspace
rm(list = ls())
require(deSolve)
# need to invert matrices
require(matlib)
library(data.table)
library(curl)

# AMY did some work here
# AMY's second change
# Amy's fourth change

# Set this working directory to where you want to save the figures.
# You need to write your own path here
###setwd('~/Desktop/MUN COVID/Second Wave/Figures/')


# Parameters as proposed for the new model

# Prob. transmission given contact (Placeholder value)
beta <- 1

# Contact rate (Placeholder value)
c <- 1

# Prob. E -> I_p given leaving E (Placeholder value)
r <- 1

# Contact tracing effectiveness rate (Placeholder value)
q <- 1

# Decrease in asymptomatic infectivity
b_a <- 0.5

# Decrease in clinical contact rate
b_c <- 0.1

# Mean time spent in Exposed
deltaE <- 1 / 4

# Mean time spent in Infected_Pre-Clinical
deltaI_p <- 1 / 2.4

# Mean time spent in Infected_Clinical
deltaI_c <- 1 / 3.2

# Mean time spent in Infected_Asymptomatic
deltaI_a <- 1 / 7

# Mean time spent in Q (Placeholder value)
deltaQ <- 1 / 2

# Mean time spent in Q_a (Placeholder value)
deltaQ_a <- 1 / 2

# Mean time spent in S_q (Placeholder value)
deltaS_q <- 1 / 2

# Days for test results (Placeholder value)
tau <- 1 / 2

# Population size (Placeholder value)
N <- 1000000


# The following is a series adjusted for quarantining given contact tracing, adjusted using Miller et. al & Tang et. al

## The following has only been partially updated (Not entirely sure options out of Q)


Miller.CT = function(t, y, parms) {
  S   <- y[1]
  E   <- y[2]
  I_p <- y[3]
  I_c <- y[4]
  I_a <- y[5]
  Q   <- y[6]
  Q_a <- y[7]
  S_q <- y[8]
  
  matrix_size <- t * 8
  
  #creating the matrix
  model_matrix <- matrix(1:matrix_size, ncol = 8)
  colnames(model_matrix) <-
    c("S", "E", "I_p", "I_c", "I_a", "Q", "Q_a", "S_q")
  
  for (i in 1:t) {
    
    model_matrix["S", i + 1] =   model_matrix["S", i]
      #got infected today
      -model_matrix["S", i] * beta * c * (model_matrix["I_p", i] + bC * model_matrix["I_c", i] + bA * model_matrix["I_a", i])
      #didn't get infected, still contacted
      -r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * (1 - beta) * q * c * (for (j in 0:tau - 1) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      } 
        b_c * x)
      + for (j in tau:(tau + 1 / deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      } 
        x)
      #finished quarantine, never had it
      +deltaS_q * model_matrix["S_q", i]
    
 
    model_matrix["E", i + 1] =   model_matrix["E", i]
      #got infected today
      +model_matrix["S", i] * beta * c * (model_matrix["I_p", i] + bC * model_matrix["I_c", i] + bA * model_matrix["I_a", i])
      #getting contacted from when they were infected
      -r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * beta * (for (j in 0:tau) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      } 
        b_c * x
      + for (j in tau:(tau + 1 / deltaE)) {
         x <- 0
         x <- x + model_matrix["S", i - j]
      }
        x)
      #going to I_p or I_a
      -deltaE * model_matrix["E", i]
    
  
    model_matrix["I_p", i + 1] =  model_matrix["I_p", i]
      #new to I_p
      + r * deltaE * model_matrix["E", i] 
      #getting contacted from when they were infected
      - r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * beta * r * (for (j in 1 / deltaE:tau + 1 / deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        b_c * x)
      #going to I_c
      - deltaI_p * model_matrix["I_p", i]
    
    
    model_matrix["I_c", i + 1] =  model_matrix["I_c", i]
      #from I_p
      + deltaI_p * model_matrix["I_p", i]
      #from Q
      + deltaQ * model_matrix["Q", i]
      #recovering
      - deltaI_c * model_matrix["I_c", i]
    
    
    model_matrix["I_a", i + 1] =  model_matrix["I_a", i]
      #new to I_a
      + (1 - r) * deltaE * model_matrix["E", i]  
      #getting contacted from when they were infected
      - r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * beta * (1 - r) * (for (j in 1 / deltaE:tau + 1 / deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        b_c * x)
      #recovering
      - deltaI_a * model_matrix["I_a", i]
    
    
    model_matrix["Q", i + 1] = model_matrix["Q", i + 1]
      #quarantined by cantact tracing
      + r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * beta * r * (for (j in 0:tau) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        b_c * x
      + for (j in tau:tau + 1/deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        x)
      #left due to (+ve) COVID test
      -deltaQ * model_matrix["Q", i]
    
      
    model_matrix["Q_a", i + 1] = model_matrix["Q_a", i + 1]
      #quarantined by cantact tracing
      + r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * beta * (1 - r) * (for (j in 0:tau) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        b_c * x
      + for (j in tau:tau + 1/deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        x)
      #finished quarantine
      -deltaQ_a * model_matrix["Q_a", i]
    
      
    model_matrix["S_q", i + 1] = model_matrix["Q_a", i + 1]
      #quarantined by cantact tracing
      + r * deltaE * model_matrix["E", i - (tau + 1 / deltaI_p)] * q * c * (1 - beta) * (for (j in 0:tau) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        b_c * x
      + for (j in tau:tau + 1/deltaI_p) {
        x <- 0
        x <- x + model_matrix["S", i - j]
      }
        x)
      #finished quarantine
      -deltaS_q * model_matrix["S_q", i]
      
  }
  
  return(list(c(dS, dE, dIP, dIC, dIA, dRS, dRA, dcumIC)))
}

# Assumed the epidemic begins with 3 exposed.
E0 = 3 / tot.popn
yini  = c(
  S = 1,
  E = E0,
  IP = 0,
  IC = 0,
  IA = 0,
  RS = 0,
  RA = 0,
  cumIC = 0
)
# the times for the numerical integration
times <- seq(0, 365, by = .1)
# The assumed R0 prior to the higher alert level
R0 <- 2.5
beta <- uniroot(find.beta, c(0, 10), R0 = R0)$root

# Day 7 for the implementation of the higher alert level
SDstart <- 7
out <- ode(
  y = yini,
  parms = NULL,
  times = times,
  func = Miller.CT
)
out7 <- data.frame(out)

# Day 10 for the implementation of the higher alert level
SDstart <- 10
out <- ode(
  y = yini,
  parms = NULL,
  times = times,
  func = Miller.CT
)
out10 <- data.frame(out)

# Making the figure for escalation on day 7 or day 10
png(file = "reescalation.png",
    width = 1000,
    height = 450)
par(mfrow = c(1, 2), mar = c(5, 6, 3, 1))
# New daily cases is the derivative of cumulative cases
plot(
  tail(out7$time, -1),
  diff(out7$cumIC) * tot.popn,
  typ = "l",
  xlab = "days",
  ylab = "daily cases",
  xlim = c(0, 30),
  main = "Re-escalation on day 7",
  ylim = c(0, 15),
  bty = "n",
  lwd = 3,
  cex.lab = 2,
  cex.main = 2,
  cex.axis = 2
)
lines(c(7, 7), c(0, max(out7$IC * tot.popn)), lty = 2, lwd = 2)
plot(
  tail(out10$time, -1),
  diff(out10$cumIC * tot.popn),
  typ = "l",
  xlab = "days",
  ylab = "daily cases",
  xlim = c(0, 30),
  main = "Re-escalation on day 10",
  ylim = c(0, 15),
  bty = "n",
  lwd = 3,
  cex.lab = 2,
  cex.main = 2,
  cex.axis = 2
)
lines(c(10, 10), c(0, max(out10$IC * tot.popn)), lty = 2, lwd = 2)
dev.off()

# Loop across escalation days
SDstart.vec <- seq(0, 20, 1)
# Preallocate output vector
out.dat = data.frame("Estart" = NULL)
Estart = NULL
ICstart = NULL
Efinal = NULL
ICfinal = NULL
LDend = NULL

# The is a loop across different days that escalation could start
for (i in seq(1, length(SDstart.vec))) {
  SDstart <- SDstart.vec[i]
  # performing the numerical integration
  out <-
    ode(
      y = yini,
      parms = NULL,
      times = times,
      func = Miller.CT
    )
  out <- data.frame(out)
  # recording the number of exposed when re-escalation starts
  Estart[i] = out$E[which(out$time == SDstart)]
  # recording the number of clinically infected when re-escalation starts
  ICstart[i] = out$IC[which(out$time == SDstart)]
  # recording the number of exposed & clinically infected at the final time
  Efinal[i] = tail(out$cumE, 1)
  ICfinal[i] = tail(out$cumIC, 1)
  # the lockdown is assumed to end when the number of clinically infected
  # is < 1. The duration of the lockdown is the end time - the start time.
  LDend[i] = out$time[max(which(out$IC * tot.popn > 1))] - SDstart
}

# Figure code
png(file = "dontwait.png",
    width = 1000,
    height = 450)
par(mfrow = c(1, 2), mar = c(5, 6, 3, 1))
plot(
  SDstart.vec,
  tot.popn * ICfinal,
  ylab = "total clinical infections",
  xlab = "Days waited till re-escalation start",
  typ = "l",
  ylim = c(0, 3000),
  bty = "n",
  lwd = 3,
  cex.lab = 2,
  cex.main = 2,
  cex.axis = 2,
  main = "Delays are bad for public health"
)
plot(
  SDstart.vec,
  LDend,
  typ = "l",
  ylab = "Days at higher alert level",
  xlab = "Days waited till re-escalation start",
  bty = "n",
  lwd = 3,
  cex.lab = 2,
  cex.main = 2,
  cex.axis = 2,
  main = "Delays are bad for well-being"
)
dev.off()
L = length(SDstart)
mod = lm(LDend[5:L] ~ SDstart.vec[5:L])
# Slope is 0.41

# Below the code investigates different R0 values
R0vec <- c(1.3, 1.2, 1.1, 1.05)
LDend = matrix(rep(0, length(SDstart.vec) * length(R0vec)),
               nrow = length(SDstart.vec),
               ncol = length(R0vec))
for (j in seq(1, length(R0vec))) {
  R0 = R0vec[j]
  beta <- uniroot(find.beta, c(0, 10), R0 = R0)$root
  for (i in seq(1, length(SDstart.vec))) {
    SDstart <- SDstart.vec[i]
    # performing the numerical integration
    out <-
      ode(
        y = yini,
        parms = NULL,
        times = times,
        func = Miller.CT
      )
    out <- data.frame(out)
    LDend[i, j] = out$time[max(which(out$IC * tot.popn > 1))] - SDstart
  }
}

png(file = "lowR0.png",
    width = 500,
    height = 450)
plot(
  SDstart.vec,
  LDend[, 1],
  typ = "l",
  ylab = "lockdown duration (days)",
  xlab = "Days waited till lockdown start",
  ylim = c(28, 38)
)
lines(SDstart.vec, LDend[, 2])
lines(SDstart.vec, LDend[, 3])
dev.off()


## Are we keeping the following? If so, how are we using it?

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

