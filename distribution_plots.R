#---------------------------------------------
  #PDF of plot (a) R code for Log-LFR distribution
  #-------------------------------------------
rm(list = ls())

x <- seq(0, 1, length.out = 100)

pva <- c(0.98, 0.85, 0.70, 0.50, 0.30)
alphava <- c(4.0, 3.0, 2.0, 0.5, 0.2)
betava <- c(1.0, 1.5, 5.0, 6.0, 7.0)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- ((p-1)*(alpha + beta*x)*W1) /
    ((1-(1-p)*W1)*(log(p)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="g(x)",
         xlim=c(0,1),
         ylim=c(0,4.5),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,1,0.2),
     labels = c("0","0.2","0.4","0.6","0.8","1"))

axis(2,
     at = seq(0.5 ,4.5,0.5),
     labels = c("0.5", "1","1.5","2","2.5","3","3.5","4","4.5"))

legend("topright",legend=c(expression(p==0.98*","~alpha==4.0*","~beta==1.0),
                              expression(p==0.85*","~alpha==3.0*","~beta==1.5),
                              expression(p==0.70*","~alpha==2.0*","~beta==5.0),
                              expression(p==0.50*","~alpha==0.5*","~beta==6.0),
                              expression(p==0.30*","~alpha==0.2*","~beta==7.0)),
              col=cols, lwd=1.5)

  #-----------------------------------------
#PDF of plot (b) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 1.5, length.out = 100)

pva <- c(0.9, 0.7, 0.5, 0.3, 0.1)
alphava <- c(0.3, 1.5, 2.0, 3.5, 4.0)
betava <- c(5.0, 4.8, 3.6, 2.4, 1.5)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- ((p-1)*(alpha + beta*x)*W1) /
    ((1-(1-p)*W1)*(log(p)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="g(x)",
         xlim=c(0,1.5),
         ylim=c(0,3),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,1.5,0.5),
     labels = c("0","0.5","1","1.5"))

axis(2,
     at = seq(0.5,3,0.5),
     labels = c("0.5", "1","1.5","2","2.5","3"))

legend("topright",legend=c(expression(p==0.9*","~alpha==0.3*","~beta==5.0),
                           expression(p==0.7*","~alpha==1.5*","~beta==4.8),
                           expression(p==0.5*","~alpha==2.0*","~beta==3.6),
                           expression(p==0.3*","~alpha==3.5*","~beta==2.4),
                           expression(p==0.1*","~alpha==4.0*","~beta==1.5)),
       col=cols, lwd=1.5)
  #-----------------------------------------
#PDF of plot (c) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 2, length.out = 100)

pva <- c(0.8, 0.6, 0.4, 0.2, 0.1)
alphava <- c(0.5, 0.5, 0.5, 0.5, 0.5)
betava <- c(2.5, 3.4, 4.6, 5.8, 7.0)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- ((p-1)*(alpha + beta*x)*W1) /
    ((1-(1-p)*W1)*(log(p)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="g(x)",
         xlim=c(0,2),
         ylim=c(0,3),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,2,0.5),
     labels = c("0","0.5","1","1.5","2"))

axis(2,
     at = seq(0.5,3,0.5),
     labels = c("0.5", "1","1.5","2","2.5","3"))

legend("topright",legend=c(expression(p==0.8*","~alpha==0.5*","~beta==2.5),
                           expression(p==0.6*","~alpha==0.5*","~beta==3.4),
                           expression(p==0.4*","~alpha==0.5*","~beta==4.6),
                           expression(p==0.2*","~alpha==0.5*","~beta==5.8),
                           expression(p==0.1*","~alpha==0.5*","~beta==7.0)),
       col=cols, lwd=1.5)
 #-----------------------------------------
#CDF plot (a) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 2, length.out = 100)
pva <- c(0.8, 0.6, 0.5, 0.4, 0.1)
alphava <- c(1, 1, 1, 1, 1)
betava <- c(0.2, 0.8, 1, 1.5, 2)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- 1-log(1-(1-p)*W1)/log(p)
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="G(x)",
         xlim=c(0,2),
         ylim=c(0,1),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,2,0.5),
     labels = c("0","0.5","1","1.5","2"))

axis(2,
     at = seq(0.1 ,1,0.1),
     labels = c("0.1", "0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))

legend("bottomright",
       legend=c("p = 0.8, α = 1, β = 0.2",
                "p = 0.6, α = 1, β = 0.8",
                "p = 0.5, α = 1, β = 1",
                "p = 0.4, α = 1, β = 1.5",
                "p = 0.1, α = 1, β = 2"),
       col=cols, lwd=1.5)
 #-----------------------------------------
#CDF plot (b) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 2, length.out = 100)
pva <- c(0.5, 0.5, 0.5, 0.5, 0.5)
alphava <- c(0.5, 1.5, 2, 2.5, 3)
betava <- c(2, 1.5, 1, 0.5, 0.2)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- 1-log(1-(1-p)*W1)/log(p)
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="G(x)",
         xlim=c(0,2),
         ylim=c(0,1),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,2,0.5),
     labels = c("0","0.5","1","1.5","2"))

axis(2,
     at = seq(0.1 ,1,0.1),
     labels = c("0.1", "0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))

legend("bottomright",
       legend=c(expression(p==0.5*","~alpha==0.5*","~beta==2),
                expression(p==0.5*","~alpha==1.5*","~beta==1.5),
                expression(p==0.5*","~alpha==2*","~beta==1),
                expression(p==0.5*","~alpha==2.5*","~beta==0.5),
                expression(p==0.5*","~alpha==3*","~beta==0.2)),
       col=cols, lwd=1.5)
  #-----------------------------------------
#HRF plot (a) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 3, length.out = 100)

pva <- c(0.9, 0.7, 0.5, 0.3, 0.1)
ava <- c(1, 1, 1, 1, 1)
bva <- c(0.01, 0.02, 0.04, 0.06, 0.08)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  a <- ava[i]
  b <- bva[i]
  
  W1 <- exp(-(a*x + (b/2)*(x^2)))
  
  h <- ((p-1)*(a + b*x)*W1) /
    ((1-(1-p)*W1)*(log(1-(1-p)*W1)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="h(x)",
         xlim=c(0,3),
         ylim=c(1,3.5),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}

axis(1,
     at = seq(0,3,0.5),
     labels = c("0","0.5","1","1.5","2","2.5","3"))

axis(2,
     at = seq(1,3.5,0.5),
     labels = c("1","1.5","2","2.5","3","3.5"))

legend("topright",
       legend=c("p = 0.9, α = 1, β = 0.02",
                "p = 0.7, α = 1, β = 0.04",
                "p = 0.5, α = 1, β = 0.06",
                "p = 0.3, α = 1, β = 0.06",
                "p = 0.1, α = 1, β = 0.08"),
       col=cols, lwd=1.5)
  #-----------------------------------------
#HRF plot (b) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 5, length.out = 100)

pva <- c(0.2, 0.4, 0.6, 0.8, 1)
alphava <- c(1.9, 1.7, 1.5, 1.3, 1.1)
betava <- c(0.1, 0.2, 0.3, 0.4, 0.5)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- ((p-1)*(alpha + beta*x)*W1) /
    ((1-(1-p)*W1)*(log(1-(1-p)*W1)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="h(x)",
         xlim=c(0,5),
         ylim=c(1,5),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}

axis(1,
     at = seq(0,5),
     labels = c("0","1","2","3","4","5"))

axis(2,
     at = seq(1,5,0.5),
     labels = c("1","1.5","2","2.5","3","3.5", "4","4.5","5"))

legend("topright",
       legend=c("p = 0.2, α = 1.9, β = 0.1",
                "p = 0.4, α = 1.7, β = 0.2",
                "p = 0.6, α = 1.5, β = 0.3",
                "p = 0.8, α = 1.3, β = 0.4",
                "p = 1, α = 1.1, β = 0.5"),
       col=cols, lwd=1.5)
  #-----------------------------------------
#HRF plot (c) R code for Log-LFR Distribution
#-------------------------------------------
  rm(list = ls())

x <- seq(0, 1, length.out = 100)

pva <- c(0.1, 0.2, 0.3, 0.4, 0.5)
alphava <- c(0.9, 0.8, 0.7, 0.6, 0.5)
betava <- c(1, 1.5, 2, 2.5, 3)

cols <- c("blue","red","cyan","magenta","green")

for(i in 1:5){
  
  p <- pva[i]
  alpha <- alphava[i]
  beta <- betava[i]
  
  W1 <- exp(-(alpha*x + (beta/2)*(x^2)))
  
  h <- ((p-1)*(alpha + beta*x)*W1) /
    ((1-(1-p)*W1)*(log(1-(1-p)*W1)))
  
  if(i==1){
    plot(x, h, type="l",
         col=cols[i], lwd=1.5,
         xlab="x", ylab="h(x)",
         xlim=c(0,1),
         ylim=c(1,4),
         xaxt="n", yaxt="n")
  } else {
    lines(x, h, col=cols[i], lwd=1.5)
  }
}
axis(1,
     at = seq(0,1,0.2),
     labels = c("0","0.2","0.4","0.6","0.8","1"))

axis(2,
     at = seq(0.5 ,4,0.5),
     labels = c("0.5", "1","1.5","2","2.5","3","3.5","4"))

legend("bottomright",
       legend=c("p = 0.1, α = 0.9, β = 1",
                "p = 0.2, α = 0.8, β = 1.5",
                "p = 0.3, α = 0.7, β = 2",
                "p = 0.4, α = 0.6, β = 2.5",
                "p = 0.5, α = 0.5, β = 3"),
       col=cols, lwd=1.5)


