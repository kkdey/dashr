

###########  Test bash - Beta Adaptive Shrinkage for smoothing  ####################


mu <- c(rep(20,1000), rep(10,1000), rep(30, 1000), rep(20, 1000))
x <- mu + rnorm(4000, 0, 1)
plot(x)

reflect <- FALSE
if(reflect){
  extended_len <- 2^{ceiling(log(length(x), base=2))}
  x_ext <- c(x, tail(x, extended_len - length(x)))
}else{
  extended_len <- 2^{ceiling(log(length(x), base=2))}
  x_ext <- c(x, rep(tail(x, 1), extended_len - length(x)))
}

pointmass=FALSE

x_odd <- x_ext[c(TRUE,FALSE)]
x_even <- x_ext[c(FALSE, TRUE)]



flow.prophat <- x_odd/(x_odd+x_even)

gridmult=sqrt(2)
mixalpha <- autoselect.mixalpha(flow.prophat = flow.prophat,
                                z_level = (x_even+x_odd),
                                mult=gridmult)
mixalpha <- mixalpha[mixalpha>0]
if(pointmass){ mixalpha = c(10^10,mixalpha) }

q <- mixalpha;

n=length(x_ext)
J=log2(n)
titable=cxxParentTItable(x_ext)
tit=titable$TItable
ptit=titable$parent

qq <- mixalpha
maxit=2000
tol=1e-5
n=dim(tit)[2]
J=dim(tit)[1]-1
nt=tit[-1,]
ns=ptit[,((1:(2*n))%%2==1)]
nf=nt-ns
maxit=maxit
tol=tol
convcrit=10.0
iter=0
lik=NULL
likk=0

pip1=0.5*rep(1,J)










z <- binshrink(x_ext, pointmass = TRUE, mode = 1)

library(smashr)
mu.est = smash(x_ext, "poiss")
