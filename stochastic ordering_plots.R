#------------------------------------------------------------------------
# Plot of an example for likelihood ratio ordering of Log-LFR distribution
#-------------------------------------------------------------------------
rm(list = ls())
alpha <- 2.5
beta  <- 1.5
y <- seq(0.01,0.99,length.out = 100)
# ---------- p1=0.15 , p2=0.2 ----------
p1 <- 0.15
A  <- (p1-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B  <- (1-((1-p1)*exp(-(alpha*log(1/(1-y))+
                         (beta/2)*(log(1/(1-y)))^2))))*log(p1)
f  <- A/B

p2 <- 0.2
A1 <- (p2-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B1 <- (1-((1-p2)*exp(-(alpha*log(1/(1-y))+
                         (beta/2)*(log(1/(1-y)))^2))))*log(p2)
g  <- A1/B1

h <- g/f

# ---------- Plot ----------
par(mgp=c(2,0.6,0))
plot(y,h,
     type="l",
     col="green",
     lwd=2,
     xlab="y",
     ylab = expression(f[Y](y) / f[X](y)),
     xaxt="n",
     yaxt="n",
     ylim=c(0.85,1.2))

# ---------- p3=0.3 , p4=0.4 ----------
p3 <- 0.3
A2 <- (p3-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B2 <- (1-(1-p3)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p3)
f1 <- A2/B2

p4 <- 0.4
A3 <- (p4-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B3 <- (1-(1-p4)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p4)
g1 <- A3/B3

h1 <- g1/f1
lines(y,h1,col="magenta",lwd=2)

# ---------- p5=0.5 , p6=0.6 ----------
p5 <- 0.5
A4 <- (p5-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B4 <- (1-(1-p5)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p5)
f2 <- A4/B4

p6 <- 0.6
A5 <- (p6-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B5 <- (1-(1-p6)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p6)
g2 <- A5/B5

h2 <- g2/f2
lines(y,h2,col="cyan",lwd=2)

# ---------- p7=0.7 , p8=0.8 ----------
p7 <- 0.7
A6 <- (p7-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B6 <- (1-(1-p7)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p7)
f3 <- A6/B6
p8 <- 0.8
A7 <- (p8-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B7 <- (1-(1-p8)*exp(-(alpha*log(1/(1-y))+
                        (beta/2)*(log(1/(1-y)))^2)))*log(p8)
g3 <- A7/B7

h3 <- g3/f3
lines(y,h3,col="black",lwd=2)
axis(1,
     at=seq(0.1,0.9,0.1))
axis(2,
     at=seq(0.9,1.1,0.05),
     las=1)
legend("bottomright",
       legend=c(
         expression(p[1]==0.1*","~~p[2]==0.2*","~~alpha==2.5*","~~beta==1.5),
         expression(p[1]==0.3*","~~p[2]==0.4*","~~alpha==2.5*","~~beta==1.5),
         expression(p[1]==0.5*","~~p[2]==0.6*","~~alpha==2.5*","~~beta==1.5),
         expression(p[1]==0.7*","~~p[2]==0.8*","~~alpha==2.5*","~~beta==1.5)
       ),
       col=c("green","magenta","cyan","black"),
       lwd=2,
       bty="o")
#-----------------------------------------------------------------------------
# Plot of counter example for likelihood ratio ordering of Log-LFR distribution
#-----------------------------------------------------------------------------
# Clear environment
rm(list=ls())

# Sequence for y
y <- seq(0.01,0.99,length.out = 100)

alpha <- 0.8
beta  <- 1.5

# -------- First Curve --------
p1 <- 0.5
A  <- (p1-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B  <- (1-((1-p1)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p1)
f  <- A/B

p2 <- 0.1
A1 <- (p2-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B1 <- (1-((1-p2)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p2)
g  <- A1/B1

h1 <- g/f
#plot
par(mgp=c(2,0.6,0))
plot(y,h1,
     type="l",
     col="green",
     lwd=2,
     xlab="y",
     ylab = expression(f[Y](y) / f[X](y)),
     xaxt="n",
     yaxt="n",
     ylim=c(0.6,2.4))

# -------- Second Curve --------
p1 <- 0.7
A  <- (p1-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B  <- (1-((1-p1)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p1)
f  <- A/B

p2 <- 0.3
A1 <- (p2-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B1 <- (1-((1-p2)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p2)
g  <- A1/B1

h2 <- g/f
lines(y, h2, col="red", lwd=2)

# -------- Third Curve --------
p1 <- 0.9
A  <- (p1-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B  <- (1-((1-p1)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p1)
f  <- A/B

p2 <- 0.6
A1 <- (p2-1)*(alpha+beta*log(1/(1-y)))*
  exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))
B1 <- (1-((1-p2)*exp(-(alpha*log(1/(1-y))+(beta/2)*(log(1/(1-y)))^2))))*log(p2)
g  <- A1/B1

h3 <- g/f
lines(y, h3, col="magenta", lwd=2)

axis(1, at=seq(0.1,0.9,0.1))
axis(2, at=seq(0.6,2.4,0.2), las=1)

legend("topright",
       legend=c(
         expression(alpha==0.8*","~~beta==1.5*","~~p[1]==0.5*","~~p[2]==0.1),
         expression(alpha==0.8*","~~beta==1.5*","~~p[1]==0.7*","~~p[2]==0.3),
         expression(alpha==0.8*","~~beta==1.5*","~~p[1]==0.9*","~~p[2]==0.6)
       ),
       col=c("green","red","magenta"),
       lwd=2,
       bty="o")
#-----------------------------------------------------------------------------
# Plot of an example for usual stochastic ordering of Log-LFR distribution
#-----------------------------------------------------------------------------
rm(list=ls())

y <- seq(0.01,0.99,length.out = 100)
p <- 0.3

# -------- 1st curve --------
alpha <- c(1,1.5,2)
beta  <- c(2.1,2.5,1.5)

alpha1 <- c(0.5,1,1.5)
beta1  <- c(2,2.5,1)

R  <- 1-log(1-(1-p)*exp(-(alpha[3]*log(1/(1-y))+(beta[3]/2)*(log(1/(1-y)))^2)))
R1 <- 1-log(1-(1-p)*exp(-(alpha1[3]*log(1/(1-y))+(beta1[3]/2)*(log(1/(1-y)))^2)))

Sol <- R1-R

plot(y,Sol,type="l",col="blue",lwd=2,
     xlab="",ylab="",
     xaxt="n",yaxt="n")   # remove default axes

# -------- 2nd curve --------
alpha <- c(2.5,3,4)
beta  <- c(2.5,3.3,4.7)

alpha2 <- c(2,2.5,3.5)
beta2  <- c(1.3,2.2,4.5)

R  <- 1-log(1-(1-p)*exp(-(alpha[3]*log(1/(1-y))+(beta[3]/2)*(log(1/(1-y)))^2)))
R1 <- 1-log(1-(1-p)*exp(-(alpha2[3]*log(1/(1-y))+(beta2[3]/2)*(log(1/(1-y)))^2)))

lines(y,R1-R,col="black",lwd=2)

# -------- 3rd curve --------
alpha <- c(5,7,9)
beta  <- c(2.7,3.8,4.3)

alpha3 <- c(4.5,6,8)
beta3  <- c(1.6,3.2,5.3)

R  <- 1-log(1-(1-p)*exp(-(alpha[3]*log(1/(1-y))+(beta[3]/2)*(log(1/(1-y)))^2)))
R1 <- 1-log(1-(1-p)*exp(-(alpha3[3]*log(1/(1-y))+(beta3[3]/2)*(log(1/(1-y)))^2)))

lines(y,R1-R,col="red",lwd=2)

# -------- 4th curve --------
alpha <- c(11,13,14.5)
beta  <- c(6.4,7.5,9.2)

alpha4 <- c(10,12,13)
beta4  <- c(4.8,5.5,8.4)

R  <- 1-log(1-(1-p)*exp(-(alpha[3]*log(1/(1-y))+(beta[3]/2)*(log(1/(1-y)))^2)))
R1 <- 1-log(1-(1-p)*exp(-(alpha4[3]*log(1/(1-y))+(beta4[3]/2)*(log(1/(1-y)))^2)))

lines(y,R1-R,col="green",lwd=2)

# -------- Axis labels --------
xlabel <-"y"
ylabel = expression(bar(F)[Y[1:n]](y)-bar(F)[X[1:n]](y))
title(xlab=" y", ylab=ylabel)
# -------- Custom axes --------
axis(1, at = seq(0.1,0.9,0.1))
axis(2, at = seq(0.00,0.10,0.02), las = 1)
# -------- Legend  --------
legend("topright",
       legend=c(
         expression(paste(alpha,"=(1,1.5,2), ",alpha^"*","=(0.5,1,1.5), ",
                          beta,"=(2.1,2.5,1.5), ",beta^"*","=(2,2.5,1)")),
         
         expression(paste(alpha,"=(2.5,3,4), ",alpha^"*","=(2,2.5,3.5), ",
                          beta,"=(2.5,4.3,6.7), ",beta^"*","=(1.3,2.2,5.8)")),
         
         expression(paste(alpha,"=(5,7,9), ",alpha^"*","=(4.5,6,8), ",
                          beta,"=(4.7,6.8,9.3), ",beta^"*","=(3.6,7.2,8.3)")),
         
         expression(paste(alpha,"=(11,13,14.5), ",alpha^"*","=(10,12,13), ",
                          beta,"=(6.4,7.5,9.2), ",beta^"*","=(4.8,5.5,8.4)"))
       ),
       col=c("blue","black","red","green"),
       lwd=2,
       bty="o",
       cex=0.8)
#----------------------------------------------------------------------------------
# Plot (a) of counter example for usual stochastic ordering of Log-LFR distribution
#----------------------------------------------------------------------------------
# Clear environment
rm(list=ls())

# Sequence for y
y <- seq(0.01,0.99,length.out=100)

# Parameter vectors
alphava <- c(0.5,1,1.5)
betava  <- c(0.1,0.3,0.5)

alpha1va <- c(2,2.5,3)
beta1va  <- c(0.2,0.4,0.6)

p <- 0.3

# Initialize vectors
R  <- rep(0,length(y))
R1 <- rep(0,length(y))

# ---- First loop ----
for(i in 1:3){
  
  alpha <- alphava[i]
  beta  <- betava[i]
  
  R <- 1 - log(1-(1-p)*exp(-(alpha*(log(1/(1-y))) +
                               (beta/2)*(log(1/(1-y)))^2)))
}

# ---- Second loop ----
for(j in 1:3){
  
  alpha1 <- alpha1va[j]
  beta1  <- beta1va[j]
  
  R1 <- 1 - log(1-(1-p)*exp(-(alpha1*(log(1/(1-y))) +
                                (beta1/2)*(log(1/(1-y)))^2)))
}

# Difference
Sol <- R1 - R

# Plot
plot(y, Sol, type="l", lty=2, col="black", lwd=2,
     xlab="y",
     ylab=expression(bar(F)[Y[1:3]](y) - bar(F)[X[1:3]](y)),
     xaxt="n")

axis(1, at=seq(0.1,0.9,0.1))


legend("topright",
       legend = expression(
        
           paste(alpha, " = (0.5, 1, 1.5), ", alpha^"*", " = (2, 2.5, 3)"),
           paste(beta, " = (0.1, 0.3, 0.5), ", beta^"*", " = (0.2, 0.4, 0.6), p = 0.3")
         
       ),
       lwd = 2,
       lty = 2,
       bty = "o")
#----------------------------------------------------------------------------------
# Plot (b) of counter example for usual stochastic ordering of Log-LFR distribution
#----------------------------------------------------------------------------------
# Clear environment
rm(list = ls())

# Sequence for y
y <- seq(0.01, 0.99, length.out = 100)

# Parameters
alphava <- c(0.1, 0.3, 0.7)
betava  <- c(2.5, 3, 3.5)

alpha1va <- c(0.6, 0.8, 1)
beta1va  <- c(1, 1.5, 2)

p <- 0.3

# ---- First loop ----
for(i in 1:3){
  alpha <- alphava[i]
  beta  <- betava[i]
  
  R <- 1 - log(1 - (1-p) * exp(-(alpha * log(1/(1-y)) +
                                   (beta/2) * (log(1/(1-y)))^2)))
}

# ---- Second loop ----
for(j in 1:3){
  alpha1 <- alpha1va[j]
  beta1  <- beta1va[j]
  
  R1 <- 1 - log(1 - (1-p) * exp(-(alpha1 * log(1/(1-y)) +
                                    (beta1/2) * (log(1/(1-y)))^2)))
}

# Difference
Sol <- R1 - R
# Plot
plot(y, Sol, type="l", lty=2, col="black", lwd=2,
     xlab="y",
     ylab=expression(bar(F)[Y[1:3]](y) - bar(F)[X[1:3]](y)),
     xaxt="n")

axis(1, at=seq(0.1,0.9,0.1))

abline(h = 0, col = "black", lwd = 2)
legend("bottomright",
       legend = expression(
         
         paste(alpha, " = (0.1, 0.3, 0.7), ", alpha^"*", " = (0.6, 0.8, 1)"),
         paste(beta, " = (2.5, 3, 3.5), ", beta^"*", " = (1, 1.5, 2), p = 0.3")
         
       ),
       lwd = 2,
       lty = 2,
       bty = "o")
#----------------------------------------------------------------------------------
# Plot (c) of counter example for usual stochastic ordering of Log-LFR distribution
#----------------------------------------------------------------------------------
# Clear environment
rm(list = ls())

# Sequence for y
y <- seq(0.01, 0.99, length.out = 100)

# Parameter vectors
alphava <- c(0.04, 0.06, 0.08)
betava  <- c(0.1, 0.3, 0.5)

alpha1va <- c(0.01, 0.02, 0.04)
beta1va  <- c(0.2, 0.4, 0.6)

p <- 0.3

# ---- First loop ----
for(i in 1:3){
  alpha <- alphava[i]
  beta  <- betava[i]
  
  R <- 1 - log(1 - (1 - p) * exp(-(alpha * log(1/(1 - y)) +
                                     (beta/2) * (log(1/(1 - y)))^2)))
}

# ---- Second loop ----
for(j in 1:3){
  alpha1 <- alpha1va[j]
  beta1  <- beta1va[j]
  
  R1 <- 1 - log(1 - (1 - p) * exp(-(alpha1 * log(1/(1 - y)) +
                                      (beta1/2) * (log(1/(1 - y)))^2)))
}

# Difference
Sol <- R1 - R

# Plot
plot(y, Sol, type="l", lty=2, col="black", lwd=2,
     xlab="y",
     ylab=expression(bar(F)[Y[1:3]](y) - bar(F)[X[1:3]](y)),
     xaxt="n")

axis(1, at=seq(0.1,0.9,0.1))

abline(h = 0, col = "black", lwd = 2)
legend("bottomleft",
       legend = expression(
         
         paste(alpha, " = (0.04, 0.06, 0.08), ", alpha^"*", " = (0.01, 0.02, 0.04)"),
         paste(beta, " = (0.1, 0.3, 0.5), ", beta^"*", " = (0.2, 0.4, 0.6), p = 0.3")
         
       ),
       lwd = 2,
       lty = 2,
       bty = "o")




