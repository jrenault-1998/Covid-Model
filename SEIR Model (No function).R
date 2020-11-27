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

#"Calculated" Beta
#beta <- 0.516

# Contact rate
#c <- 5

# Prob. transmission given contact
#alpha = beta/c
alpha.vec <- seq(0.1,0.2,0.1)

# Days the model runs
t <- 100

# Days for test results
tau <- 2

# Prob. E -> I_p given leaving E
r <- mean(0.4,0.25,0.37,0.42,0.51,0.59,0.2,0.76,0.76)

# Contact tracing effectiveness rate
q.vec <-  seq(0,1,0.5) #0.75

# Decrease in asymptomatic infectivity
b_a <- 0.5

# Decrease in clinical contact rate
b_c <- 0.1

# Mean time spent in Exposed
deltaE <- 1 / 3

# Mean time spent in Infected_Pre-Clinical
deltaI_p <- 1 / 2

# Mean time spent in Infected_Clinical
deltaI_c <- 1 / 3

# Mean time spent in Infected_Asymptomatic
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


timeloop = function(c){
y  = c(
  S = (N - E0),
  E = E0,
  I_p = 0,
  I_c = 0,
  I_a = 0,
  Q = 0,
  Q_a = 0,
  S_q = 0,
  R = 0
)


CT_size <- (t+1) * 8

#creating the matrix
CT <- matrix(0, nrow = t+1, ncol = 9)
colnames(CT) <- c("S", "E", "I_p", "I_c", "I_a", "Q", "Q_a", "S_q", "R")
CT[1, "S"] = N
CT[2, "S"] = N
CT[3, "S"] = N
CT[4, "S"] = N
CT[5, "S"]   = y[1]
CT[5, "E"]   = y[2]
CT[5, "I_p"] = y[3]
CT[5, "I_c"] = y[4]
CT[5, "I_a"] = y[5]
CT[5, "Q"]   = y[6]
CT[5, "Q_a"] = y[7]
CT[5, "S_q"] = y[8]
CT[5, "R"]   = y[9]
CT[1, "E"]   = 0
CT[2, "E"]   = 0
CT[3, "E"]   = 0
CT[4, "E"]   = 0

output <- matrix(0,nrow = length(q.vec), ncol=length(alpha.vec))

# for(k in 1:length(alpha.vec)){
#   alpha<-alpha.vec[k]
#   
#   for(j in 1:length(q.vec)){
#     
#     q<-q.vec[j]
    
    for (i in seq(5,t-1) ) {
      S   = CT[i,"S"]
      S1  = CT[i-1,"S"]
      S2  = CT[i-2,"S"]
      S3  = CT[i-3,"S"]
      E   = CT[i,"E"]
      E4  = CT[i-4,"E"]
      I_p = CT[i,"I_p"]
      I_c = CT[i,"I_c"]
      I_a = CT[i,"I_a"]
      Q   = CT[i,"Q"]
      Q_a = CT[i,"Q_a"]
      S_q = CT[i,"S_q"]
      R   = CT[i, "R"]
      N   = S + E + I_p + I_c + I_a + Q + Q_a + S_q + R
      
      CT[i+1, "S"] =   S - S * alpha * c * (I_p + b_c * I_c + b_a * I_a)/N- (E4*r*deltaE) * q * (1-alpha) * c * (b_c*(S+S1)+(S2+S3))/N + deltaS_q * S_q
      
      
      CT[i+1, "E"] =   E + S * alpha * c * (I_p + b_c * I_c + b_a * I_a)/N- (E4*r*deltaE) * c * q * alpha * (b_c*(S+S1)+(S2+S3))/N- deltaE * E
      
      
      CT[i+1, "I_p"] =  I_p+ r * deltaE * E - deltaI_p * I_p
      
      
      CT[i+1, "I_c"] =  I_c+ deltaI_p * I_p+ deltaQ * Q- deltaI_c * I_c
      
      
      CT[i+1, "I_a"] =  I_a+ (1-r) * deltaE * E- deltaI_a * I_a
      
      
      CT[i+1, "Q"] =  Q+ (E4*r*deltaE) * c * q * alpha * r * (b_c*(S+S1)+(S2+S3))/N- deltaQ * Q
      
      
      CT[i+1, "Q_a"] =  Q_a+ (E4*r*deltaE) * c * q * alpha * (1-r) * (b_c*(S+S1)+(S2+S3))/N-deltaQ_a * Q_a
      
      
      CT[i+1, "S_q"] =  S_q + (E4*r*deltaE) * q * (1-alpha) * c * (b_c*(S+S1)+(S2+S3))/N - deltaS_q * S_q
      
      CT[i+1, "R"] =  R + deltaI_c * I_c + deltaQ_a * Q_a + deltaI_a * I_a
    }
    
diff.target = unname(CT[i+1,"R"]) - target
# The function returns the differnce between the final value of R and the target value
return(diff.target)
    
#  }}    #Commented for loop
}

#preallocate the output
output = matrix(rep(0,length(alpha.vec)*length(q.vec)),nrow=length(q.vec), ncol=length(alpha.vec))
# specify start and end values for the search across c values.
c.start <- 0
c.end <- 20

for(k in 1:length(alpha.vec)){
  alpha<-alpha.vec[k]
  for(j in 1:length(q.vec)){
    q <-q.vec[j]
    target <- 0.001
     # The if clause is to prevent uniroot giving an error if there
     # exists no value of the final R such that the target is meet.
     # if(sign(timeloop(c.start))==sign(timeloop(c.end))){
     #   a <- 1 #output[j,k]=NA
     #   }
     #else{
    print(c(alpha,q))
      cval <- uniroot(timeloop,c(c.start,c.end))$root
      output[j,k] <- cval
    
      #}
  }}

filled.contour(q.vec,alpha.vec,output, xlab = "q", ylab = "alpha")




# plot (q.vec, output[,1], typ = "l", ylim = c(0, 1), xlim = c(0, 1), ylab = "R Final", xlab = "Contact Tracing Effectiveness", main = "c=5")
# lines(q.vec, output[,2], col = "red")
# lines(q.vec, output[,3], col = "blue")
# lines(q.vec, output[,4], col = "orange")
# lines(q.vec, output[,5], col = "green")
# legend( "topright", c("alpha=0.12", "alpha=0.14", "alpha=0.16", "alpha=0.18", "alpha=0.20"),
#                 text.col=c("black", "red", "blue", "orange", "green") )

# df <- data.frame("S" = CT[, "S"], "E" = CT[, "E"],"I_p" = CT[, "I_p"], "I_c" = CT[,"I_c"],
#                  "I_a" = CT[,"I_a"], "Q" = CT[, "Q"], "Q_a" = CT[,"Q_a"], "S_q" = CT[,"S_q"], "R" = CT[, "R"])
# 
