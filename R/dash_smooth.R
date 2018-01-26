#' @title Adaptive smoothing using Beta adaptive shrinkage
#'
#' @description Given a vector of counts, which are noisy estimates of an
#' underlying Poisson counts data, the function performs adaptive smoothing
#' of the counts by fitting a Beta adaptive shrinkage model.
#'
#' @details The input to dash-smooth is a vector of counts which are noisy
#' versions of a smooth process. We fit a multiscale model on these counts
#' and the message flow proportions are assumed to have a flexible mixture
#' Beta prior centered around mean of 0.5.
#'
#' @param x, a vector of counts
#' @param concentration a vector of concentration scales for different Dirichlet
#'                   compositions. Defaults to NULL, in which case, we append
#'                   concentration values of Inf, 100, 50, 20, 10, 5, 2, 1, 0.5
#'                   and 0.1.
#' @param pi_init An initial starting value for the mixture proportions. Defaults
#'                to same proportion for all categories.
#' @param reflect Boolean indicating if the vector is padded by a reflection of the tail
#'                or the tailmost value so that the padded vector has length a power of 2.
#' @param squarem_control A list of control parameters for the SQUAREM/IP algorithm,
#'              default value is set to be control.default=list(K = 1, method=3,
#'               square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1,
#'               objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE).
#' @param dash_control A list of control parameters for determining the concentrations
#'                     and prior weights and fdr control parameters for dash fucntion.
#' @return Returns a list of the following items
#'                     \code{estimate}: The adaptively smoothed values of the counts vector \code{x}.
#'                     \code{pi_weights}: The mixture proportions estimated from different levels of multiscale model.
#'                     \code{loglik}: The loglikelihood value of the fitted model.
#'
#' @examples
#' mu <- c(rep(10, 50), rep(20, 50), rep(30, 50), rep(10, 50))
#' x <- mu + rnorm(200, 0, 1)
#' out <- dash_smooth(x)
#' out$estimate
#'
#' @export


dash_smooth = function(x, concentration = NULL,
                       pi_init = NULL, reflect = FALSE,
                       squarem_control = list(),
                       dash_control = list()){

  require(Rcpp)
  require(inline)

  if(min(x) < 0){stop ("negative values in x not permitted")}

  dash_control.default <- list(add_NULL = TRUE, add_Inf = TRUE, add_corner = FALSE,
                               corner_val = 0.005, null_weight = 1, Inf_weight = 100,
                               corner_weight = 1, fdr_bound = 50)

  dash_control <- modifyList(dash_control.default, dash_control)

  squarem_control.default=list(K = 1, method=3,
                               square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1,
                               objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)

  squarem_control <- modifyList(squarem_control.default, squarem_control)


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
  concentration <- unique(concentration)


  conc <- concentration
  conc[conc == Inf] <- 10^10

  mode <- rep(1, 2)
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

  if(reflect){
    extended_len <- 2^{ceiling(log(length(x), base=2))}
    x_ext <- c(x, tail(x, extended_len - length(x)))
  }else{
    extended_len <- 2^{ceiling(log(length(x), base=2))}
    x_ext <- c(x, rep(tail(x, 1), extended_len - length(x)))
  }

  x_odd <- x_ext[c(TRUE,FALSE)]
  x_even <- x_ext[c(FALSE, TRUE)]

  titable=cxxParentTItable(x_ext)

  tit=titable$TItable
  ptit=titable$parent
  n=dim(tit)[2]
  J=dim(tit)[1]-1
  nt=tit[-1,]
  ns=ptit[,((1:(2*n))%%2==1)]
  nf=nt-ns
  loglik = 0

  pi_weights <- matrix(0, dim(ns)[1], length(concentration))
  pb <- txtProgressBar(min = 0, max = dim(ns)[1], style = 3)

  for(k in 1:dim(ns)[1]){
    mat <- rbind(ns[k,], nf[k, ])
    mat <- mat + 1
    out <- dash2(mat, conc_mat, prior, pi_init, squarem_control)
    loglik = loglik + out$loglik
    pi_weights[k,] <- out$pi
    setTxtProgressBar(pb, k)
  }

  close(pb)

  concentration2 <- concentration
  concentration2[concentration2 == Inf] <- 10^10
  est=reverse.pp(tit,pi_weights,concentration2,2)
  est2 <- est - mean(est) + mean(x_ext)

  ll <- list("estimate" = est2[1:length(x)],
               "pi_weights" = pi_weights,
               "loglik" = loglik)

  return(ll)
}


###############################################################################
###############################################################################
###############a skeleton version of the dash() function  ######################
###############################################################################
###############################################################################



dash2 <- function(comp_data,
                  conc_mat,
                  prior,
                  pi_init,
                  squarem_control){
  comp_data <- t(comp_data)
  matrix_log_lik <- matrix(0, dim(comp_data)[1], dim(conc_mat)[1])

  for(n in 1:dim(comp_data)[1]){
    x <- comp_data[n,]
    for(k in 2:dim(conc_mat)[1]){
      # numero <- sum(x)*beta(sum(conc_mat[k,]), sum(x))
      lognumero <- log(sum(x)) - LaplacesDemon::ddirichlet(rep(1,2), alpha = c(sum(conc_mat[k,]), sum(x)), log=TRUE)
      if(lognumero == -Inf | lognumero == Inf ){
        matrix_log_lik[n, k] <- lognumero
      }else{
        index1 <- which(x > 0)
        logdeno <- sum(log(x[index1]) -  sapply(1:length(index1), function(mm) return(LaplacesDemon::ddirichlet(rep(1,2), alpha = c(conc_mat[k, index1[mm]], x[index1[mm]]), log=TRUE))))
        matrix_log_lik[n,k] <- lognumero - logdeno
      }
    }
    matrix_log_lik[n,1] <- logfac(sum(x)) - sum(sapply(x, function(y) return(logfac(y)))) + sum(x*log((conc_mat[1,]+1e-04)/sum(conc_mat[1,]+1e-04)))
  }

  matrix_lik <- exp(matrix_log_lik - apply(matrix_log_lik, 1, function(x) return(max(x))) %*% t(rep(1, dim(matrix_log_lik)[2])))
  fit=do.call("mixEM",args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=squarem_control))
  pihat <- fit$pihat

  lll <- list ("pi" = pihat,
               "loglik" = sum(matrix_log_lik))
  return(lll)
}


interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


#############  Shift operator functions  ##########################

rshift = function(x){L=length(x); return(c(x[L],x[-L]))}

lshift = function(x){return(c(x[-1],x[1]))}


##############  Parent TI table maker (R version)  ######################

ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)

  # Create decomposition table of signal, using pairwise sums,
  # keeping just the values that are *not* redundant under the
  # shift-invariant scheme.  This is very similar to TI-tables
  # in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table

  for(D in 0:(J-1)){
    nD = 2^(J-D);
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}

#################  Parent TI table (C++ version)  #######################

src <- '
        NumericVector signal=sig;
        int n=(int) signal.size();
        int J=(int) log2((double)n);

        NumericMatrix parent(J,2*n);
        NumericMatrix TItable(J+1,n);
        TItable(0,_) = signal;
        for (int D=0; D<J; D++){
        int nD=(int) pow(2., (double) (J-D)), pD=(int) pow(2.,(double) D);
        for (int l=0; l<pD; l++){
        int a=l*nD+1, b=2*l*nD+1, d;
        for (int i=0; i<nD-1; i++){
        d=TItable(D,a+i-1);
        parent(D,b+i-1)=d;
        parent(D,b+i+nD)=d;
        }
        //i=nD-1
        d=TItable(D,a+nD-2);
        parent(D,b+nD-2)=d;
        parent(D,b+nD-1)=d;

        for (int i=0; i<nD; i++)
        TItable(D+1,a+i-1)=parent(D,b+2*i-1)+parent(D,b+2*i);
        }
        }
        return(List::create(Named("TItable")=TItable, Named("parent")=parent));
'
cxxParentTItable <- cxxfunction(signature(sig="numeric"),
                                body=src,
                                plugin="Rcpp",
                                inc="#include <cmath>")


#This function computes the posterior means
sfunc=function(p,q,x0,x1,mode){
  if(mode==1){
    xx = x0 + x1
    if(p==1){
      ss = 0.5*rep(1,length(x0))
    }else if(p == 0){
      ss = (x0 + q)/(xx + 2*q)
    }else{
      # Compute first half denomenator of sum.
      fhd = log(1-p) - log(p) + (xx+1)*log(2) + lbeta(x0+q,x1+q) - lbeta(q,q)
      fhd = 2 + exp(fhd)

      # Compute second half denomenator.
      shd1 = log(p) - log(1-p) + lbeta(q,q) - xx*log(2) - lbeta(x0+q+1,x1+q)
      shd  = exp(shd1) + (xx + 2*q)/(x0 + q)

      # Put together two pieces.  Numerators are just 1.
      ss = (1/fhd) + (1/shd)
    }
  }else if(mode==2){
    nq = length(q)
    nn = length(x0)
    x0m = rep(1,nq)%o%x0
    x1m = rep(1,nq)%o%x1
    pm = p%o%rep(1,nn)
    qm = q%o%rep(1,nn)
    rq = rep(1,nq)
    num = log(pm)+lbeta(x0m+qm+1,x1m+qm)-lbeta(qm,qm)
    numpm = apply(num,2,max)
    num = num-rq%o%numpm
    num = log(colSums(exp(num)))+numpm
    den = log(pm)+lbeta(x0m+qm,x1m+qm)-lbeta(qm,qm)
    denpm = apply(den,2,max)
    den = den-rq%o%denpm
    den = log(colSums(exp(den)))+denpm
    ss = exp(num-den)
  }
  return(ss)
}


################  Getting dash shrunk estimates   ##########################


#This is similar to reverse.pwave used in cyclespin.smooth
reverse.pp=function(dmat,pp,qq,mode){
  n=dim(dmat)[2]
  J=log2(n)
  if(mode==1){
    qq=rep(qq,J)
  }else if(mode==2){
    qq=rep(1,J)%o%qq
  }
  # Beginning with the total number of counts, working from coarse
  # scales (i.e., deep depths) upwards, gradually build up the MSPB
  # shrinkage estimate by multiplying by appropriate shrinkage factors
  # at each level.
  est = dmat[J+1,]
  for (D in J:1){
    nD = 2^(J-D+1)
    nDo2 = nD/2
    for (l in 0:(2^(D-1)-1)){
      # Set indexing so as to pick off blocks of size 2^(J-D+1)
      # when shrinking estimates at depth D+1 down to finer
      # scale at depth D.
      ind = (l*nD+1):((l+1)*nD)
      xxD = dmat[D,ind]
      estvec = est[ind]

      # In the first half of the vector of D+1-depth estimates,
      # we can shrink using the D-depth counts in the order
      # in which they appear.
      estl = estvec[1:nDo2]
      xxDl = xxD[seq(1,nD-1,2)]
      xxDr = xxD[seq(2,nD,2)]
      if(mode==1){
        ss = sfunc(pp[D],qq[D],xxDl,xxDr,1)
      }else if(mode==2){
        ss = sfunc(pp[D,],qq[D,],xxDl,xxDr,2)
      }
      nestl = interleave(estl*ss,estl*(1 - ss))

      # In the second half of the vector of D+1-depth counts,
      # we right-shift the D-depth counts, compute the shrunken
      # values, and left-shift these values back to the order
      # of the above.
      estr = estvec[(nDo2+1):nD]
      sxxD = rshift(xxD)
      xxDl = sxxD[seq(1,nD-1,2)]
      xxDr = sxxD[seq(2,nD,2)]
      if(mode==1){
        ss = sfunc(pp[D],qq[D],xxDl,xxDr,1)
      }else if(mode==2){
        ss = sfunc(pp[D,],qq[D,],xxDl,xxDr,2)
      }
      nestr = interleave(estr*ss,estr*(1 - ss))
      nestr = lshift(nestr)

      # Combine the estimates from both halves of the D+1-depth
      # counts, and store.
      est[ind] = 0.5*( nestl + nestr )
    }
  }
  return(est)
}
