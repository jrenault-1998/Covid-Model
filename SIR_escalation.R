# Final size
R1 = 3
R2 = 0.4
gamma = 1/13
I0 = 1

# Time at peak
te = seq(0.1,30,.1)
# Active infecteds at peak
I.te = I0*exp(gamma*(R1-1)*te)
# Cumulative infected at peak
C.minus = R1*I0*(exp(gamma*(R1-1)*te)-1)/(R1-1)
# Cumulative new infections after escalation
C.plus = R2*I.te/(1-R2)

## Example figure
te.val = 10
# Check height at peak
Ival = I0*exp(gamma*(R1-1)*te.val)
tvec1 = seq(0,te.val,.1)
tvec2 = seq(te.val,100,.1)
plot(tvec1, I0*exp((gamma*(R1-1)*tvec1)), col = "red", typ="l",xlim = c(0,100), ylab ="Infected", ylim=c(0,Ival), xlab = "day")
lines(tvec2, I0*exp(gamma*(R1-1)*te.val)*exp(gamma*(R2-1)*(tvec2-te.val)), col = "blue")

# Time vs. final size
plot(te, C.minus, typ = "l", col = "red", ylab = "cumulative cases before and after peak", xlab = "escalation day")
lines(te, C.plus, col = "blue")

# Cumulative infections vs. escalation day
plot(te, C.minus+C.plus, typ="l", ylab = "final size", xlab = "cumulative inf on escalation day")

# Cumulative infections at peak vs. final size
plot(C.minus, C.minus+C.plus, typ="l", ylab = "final size", xlab = "cumulative inf on escalation day")

final.size=C.minus+C.plus
mod=lm(final.size~C.minus)
print(summary(mod))

predicted.slope = 1 + R2*(R1-1)/(1-R2)/R1
predicted.intercept = R2/(1-R2)
