# clear workspace
rm(list = ls())
require(deSolve)

# Set this working directory to where you want to save the figures.
# You need to write your own path here
###setwd('~/Desktop/MUN COVID/Second Wave/Figures/')

# Parameters as proposed for the new model

#"Calculated" Beta
#beta <- 0.516

## AH: Commented out, c is now something to solve for.
# Contact rate
#c.vec <- seq(2,4,0.5)

# Prob. transmission given contact
#alpha = beta/c
#alpha <- 0.18
alpha.vec <- seq(0.1,0.2,0.02)

# Days the model runs
t <- 100

# Days for test results
tau <- 2

# Prob. E -> I_p given leaving E #Published Miller (r=0.4)
r <- mean(0.4,0.25,0.37,0.42,0.51,0.59,0.2,0.76,0.76)

# Contact tracing effectiveness rate
q.vec <-  seq(0,1,0.2) #0.75

### AH: Please check these numbers with the published Miller,
### Metcalfe, Grenfell paper published in NAture Medicine - the
### numbers changed in the preprint vs. the published version.

# Decrease in asymptomatic infectivity
b_a <- 0.5

# Decrease in clinical contact rate      ##Published Miller(consider changing to vec of values? (0.1,0.5,1))
b_c <- 0.1

# Mean time spent in Exposed             #Published Miller
deltaE <- 1 / 3

# Mean time spent in Infected_Pre-Clinical
deltaI_p <- 1 / 2

# Mean time spent in Infected_Clinical
deltaI_c <- 1 / 3

# Mean time spent in Infected_Asymptomatic   #Published Miller
deltaI_a <- 1 / 5

# Mean time spent in Q
deltaQ <- 1 / ((1/deltaE + 1/deltaI_p)/2 + tau)

# Mean time spent in Q_a
deltaQ_a <- 1 / 14

# Mean time spent in S_q
deltaS_q <- 1 / 14

# Initial Population size
N <- 1

# Assumed the epidemic begins with 0.0001% of population exposed
# Canada's March 11 confirmed cases of covid
E0 <- 2.74*10^(-6)

# Returns the value of R at the final time for a given value of c.
timeloop = function(c){
  # Initial values:
  S = (N - E0)
  E = E0
  I_p = 0
  I_c = 0
  I_a = 0
  Q = 0
  Q_a = 0
  S_q = 0
  R = 0
  # Enter the values into df for t = 0 to 4 - establishes values for first 5 days
  df <- data.frame(time = seq(0:4), S=S, E=E, I_p=I_p, I_c = I_c, I_a = I_a, Q=Q, Q_a = Q_a, S_q = S_q, R=R)
  for (i in seq(5,t-1) ) {
    # Lagged variables
    S1  = df[i-1,"S"]
    S2  = df[i-2,"S"]
    S3  = df[i-3,"S"]
    E   = df[i,"E"]
    E4  = df[i-4,"E"]
    I_p = df[i,"I_p"]
    I_c = df[i,"I_c"]
    I_a = df[i,"I_a"]
    Q   = df[i,"Q"]
    Q_a = df[i,"Q_a"]
    S_q = df[i,"S_q"]
    R   = df[i, "R"]
    N   = S + E + I_p + I_c + I_a + Q + Q_a + S_q + R
    
    df[i+1, "S"] =   S - S * alpha * c * (I_p + b_c * I_c + b_a * I_a)/N- (E4*r*deltaE) * q * (1-alpha) * c * (b_c*(S+S1)+(S2+S3))/N + deltaS_q * S_q
    
    
    df[i+1, "E"] =   E + S * alpha * c * (I_p + b_c * I_c + b_a * I_a)/N- (E4*r*deltaE) * c * q * alpha * (b_c*(S+S1)+(S2+S3))/N- deltaE * E
    
    
    df[i+1, "I_p"] =  I_p+ r * deltaE * E - deltaI_p * I_p
    
    
    df[i+1, "I_c"] =  I_c+ deltaI_p * I_p+ deltaQ * Q- deltaI_c * I_c
    
    
    df[i+1, "I_a"] =  I_a+ (1-r) * deltaE * E- deltaI_a * I_a
    
    
    df[i+1, "Q"] =  Q+ (E4*r*deltaE) * c * q * alpha * r * (b_c*(S+S1)+(S2+S3))/N- deltaQ * Q
    
    
    df[i+1, "Q_a"] =  Q_a+ (E4*r*deltaE) * c * q * alpha * (1-r) * (b_c*(S+S1)+(S2+S3))/N-deltaQ_a * Q_a
    
    
    df[i+1, "S_q"] =  S_q + (E4*r*deltaE) * q * (1-alpha) * c * (b_c*(S+S1)+(S2+S3))/N - deltaS_q * S_q
    
    df[i+1, "R"] =  R + deltaI_c * I_c + deltaQ_a * Q_a + deltaI_a * I_a
  }
  diff.target = df[i+1,"R"] - target
  # The function returns the differnce between the final value of R and the target value
  return(diff.target)
}

# The following is a series adjusted for quarantining given contact tracing, adjusted using Miller et. al & Tang et. al

## The following has only been partially updated (Not entirely sure options out of Q)

#preallocate the output
output = matrix(rep(0,length(alpha.vec)*length(q.vec)),nrow=length(q.vec), ncol=length(alpha.vec))
# specify start and end values for the search across c values.
c.start <- 0
c.end <- 20
target <- 0.001

for(k in 1:length(alpha.vec)){
  alpha<-alpha.vec[k]
  for(j in 1:length(q.vec)){
    q<-q.vec[j]
    # The if clause is to prevent uniroot giving an error if there
    # exists no value of the final R such that the target is meet.
    if(sign(timeloop(c.start))==sign(timeloop(c.end))){
      output[j,k]=NA}
    else{
      cval <- uniroot(timeloop,c(c.start,c.end))$root
      output[j,k] <- cval}
  }}

filled.contour(q.vec,alpha.vec,output, xlab = "q", ylab = "alpha")


for(k in 1:length(alpha.vec)){
  alpha<-alpha.vec[k]
  for(j in 1:length(q.vec)){
    q<-q.vec[j]

# However, I am not clear why R final can be bigger than 1. I think there
# is a mistake in your original code, which I just put inside a function
# and left unchanged. If N = 1 and there are no births and deaths, then
# R cannot be bigger than 1, because the sum of all the variables has to be 1 for all time.
