
########### test dash_smooth()  using the core function #################


mu <- c(rep(10, 1000), rep(20, 1000), rep(30, 1000), rep(10, 1000))
x <- mu + rnorm(4000, 0, 1)

##########       dashu         #############################

out <- dashu(x)

##########    smash   ############################

smash_out <- smashr::smash.poiss(x)

################  Visualization  ######################

plot(x, col = "gray80")
lines(mu, col = "blue", lwd = 2)
lines(out$estimate, col = "red", lwd = 2)
lines(smash_out, col = "yellow")
legend("topright", # places a legend at the appropriate place
       c("truth","dashu", "smash"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       cex = 0.6,
       col=c("blue","red", "yellow"))



mu <- c(rep(10, 1000), rep(20, 1000), rep(30, 1000), rep(10, 1000))
x <- mu + c(rnorm(1000, 0, 0.01), rnorm(1000, 0, 1), rnorm(1000, 0, 5), rnorm(1000, 0, 10))
x[x < 0] = 0
y <- sapply(x, function(z) return(rpois(1, z)))

##########       dashu         #############################

out <- dashu(y)

##########    smash   ############################

smash_out <- smashr::smash.poiss(y)


################  Visualization  ######################

plot(y, col = "gray80")
lines(out$estimate, col = "red", lwd = 2)
lines(smash_out, col = "yellow")
lines(mu, col = "blue", lwd = 2)
legend("topright", # places a legend at the appropriate place
       c("truth","dashu", "smash"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       cex = 0.6,
       col=c("blue","red", "yellow"))




mu <- c(rep(10, 1000), rep(20, 1000), rep(30, 1000), rep(10, 1000))
x <- mu + c(rnorm(1000, 0, 1), rnorm(1000, 0, 1), rnorm(1000, 0, 1), rnorm(1000, 0, 1))
x[x < 0] = 0
y <- sapply(x, function(z) return(rpois(1, z)))

##########       dashu         #############################

out <- dashu(y)

##########    smash   ############################

smash_out <- smashr::smash.poiss(y)


################  Visualization  ######################

plot(y, col = "gray80")
lines(smash_out, col = "yellow", lwd = 2)
lines(out$estimate, col = "red", lwd = 2)
lines(mu, col = "blue", lwd = 2)
legend("topright", # places a legend at the appropriate place
       c("truth","dashu", "smash"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       cex = 0.6,
       col=c("blue","red", "yellow"))



mu <- c(rep(10, 1000), rep(20, 1000), rep(30, 1000), rep(10, 1000))
x <- mu + c(rnorm(1000, 0, 10), rnorm(1000, 0, 10), rnorm(1000, 0, 10), rnorm(1000, 0, 10))
x[x < 0] = 0
y <- sapply(x, function(z) return(rpois(1, z)))

##########       dashu         #############################

out <- dashu(y)

##########    smash   ############################

smash_out <- smashr::smash.poiss(y)


################  Visualization  ######################

plot(y, col = "gray80")
lines(smash_out, col = "yellow", lwd = 2)
lines(out$estimate, col = "red", lwd = 2)
lines(mu, col = "blue", lwd = 2)
legend("topright", # places a legend at the appropriate place
       c("truth","dashu", "smash"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       cex = 0.6,
       col=c("blue","red", "yellow"))

