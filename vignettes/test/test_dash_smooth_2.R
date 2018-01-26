

mu <- c(rep(10, 1000), rep(20, 1000), rep(30, 1000), rep(10, 1000))
x <- mu + rnorm(4000, 0, 1)


input <- get(load("../test/reads_all_1_93297593_93307481.Robj"))
x <- input[[1]][[2]][1,]
reflect <- FALSE
ll <- length(x)
if(reflect){
  extended_len <- 2^{ceiling(log(length(x), base=2))}
  x_ext <- c(x, tail(x, extended_len - length(x)))
}else{
  extended_len <- 2^{ceiling(log(length(x), base=2))}
  x_ext <- c(x, rep(tail(x, 1), extended_len - length(x)))
}


x_odd <- x_ext[c(TRUE,FALSE)]
x_even <- x_ext[c(FALSE, TRUE)]

######################### create TI tables for the observations  ########################

titable=cxxParentTItable(x_ext)
tit=titable$TItable
ptit=titable$parent
n=dim(tit)[2]
J=dim(tit)[1]-1
nt=tit[-1,]
ns=ptit[,((1:(2*n))%%2==1)]
nf=nt-ns

dash_control.default <- list(add_NULL = TRUE, add_Inf = TRUE, add_corner = FALSE,
                             corner_val = 0.005, null_weight = 1, Inf_weight = 1,
                             corner_weight = 1, fdr_bound = 50)

dash_control <- modifyList(dash_control.default, dash_control)

squarem_control.default=list(K = 1, method=3,
                             square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1,
                             objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)

squarem_control <- modifyList(squarem_control.default, squarem_control)

## check if the mode has same length as columns of composition data

if(verbose){
  cat("Checking inputs and processing the data \n")
}

############  Determine the mode of the Dirichlet distribution   ##################

if(is.null(mode)){
  mode <- rep(1, dim(comp_data)[2])
}else{
  mode <- mode/min(mode[mode!=0])
}

## weights for the samples - check if this vector has the same
## length as the number of samples (number of rows of the compositional data),
## unless it is NULL, in which case, all samples have equal weights.

if(!is.null(sample_weights)){
  if(length(sample_weights) != dim(comp_data)[1]){
    stop("The length of the user-defined sample weights must match with number of rows
         in the comp_data")
  }
  }

## check if an initial mixture proportion pi has been provided by the user.

if(!is.null(pi_init)){
  if(length(pi_init) != dim(comp_data)[2]){
    stop("The length of the user-defined pi_init must match with number of columns
         in the comp_data")
  }
  }

## if background mode or probability is provided, we check if the length of the
## vector matches with number of samples

if(!is.null(mode)){
  if(length(mode) != dim(comp_data)[2]){
    stop("The length of the user-defined mode must match with number of columns
         in the comp_data")
  }
  }

## add prior concentrations - adding an Inf and 1 concentration to the mix
## if not provided by the user

concentration <- unique(concentration)

if(is.null(concentration)){
  concentration <- c(Inf, 100, 50, 20, 10, 5, 2, 1)
}else{
  if (dash_control$add_NULL){
    concentration <- c(concentration, 1)
  }
  if (dash_control$add_Inf){
    concentration <- c(concentration, Inf)
  }
  if(dash_control$add_corner){
    if (min(concentration) > dash_control$corner_val){
      concentration <- c(concentration, dash_control$corner_val)
    }
  }
  concentration <- sort(concentration, decreasing = TRUE)
}

conc <- concentration
conc[conc == Inf] <- 10^5
conc_mat <- t(sapply(conc,function(x) return(x*(mode+1e-04))))

prior <- array(1, length(concentration))
if(length(which(concentration == Inf)) > 0){
  prior[which(concentration == Inf)] <- dash_control$Inf_weight
}
if(length(which(concentration == 1)) > 0){
  prior[which(concentration == 1)] <- dash_control$null_weight
}
if(min(concentration) < 1){
  prior[which(concentration == min(concentration))] <- dash_control$corner_weight
}



pi_weights <- matrix(0, dim(ns)[1], length(concentration))
for(k in 1:dim(ns)[1]){
  mat <- rbind(ns[k,], nf[k, ])
  mat <- mat + 1
  out <- dash(mat, optmethod = "mixEM", verbose=TRUE)
  pi_weights[k,] <- out$fitted_pi
}

concentration2 <- concentration
concentration2[concentration2 == Inf] <- 10^10
est=reverse.pp(tit,pi_weights,concentration2,mode)

est_mod <- est + mean(x_ext) - mean(est)
plot(est_mod, type = "l", col = "blue")
lines(x_ext, col = "red")
lines(est3, col = "green")
est3 <- smashr::smash.poiss(x_ext)




mat <- rbind(c(5, 0, 2, 0),
             c(1, 1, 0, 1),
             c(100, 100, 50, 100),
             c(20, 50, 100, 10),
             c(10, 10, 200, 20),
             c(50, 54, 58, 53),
             c(1,1,1,3),
             c(2, 4, 1, 1))
out <- dash(mat, optmethod = "mixEM", verbose=TRUE)
out <- dash(mat, optmethod = "w_mixEM", verbose=TRUE)








