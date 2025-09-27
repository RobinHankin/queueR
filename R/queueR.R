#' @useDynLib queueR, .registration=TRUE

#' @importFrom stats integrate
#' @importFrom stats qgamma
#' @importFrom Rcpp evalCpp


#' @export
poch_complex_n <- function(a,q=a,alpha,maxit=1e6){
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  numerator <- jj[,1] + 1i*jj[,2]

  a <- a*q^alpha
  jj <- cbind(Re(a),Im(a),Re(q),Im(q),Inf)
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  denominator <- jj[,1] + 1i*jj[,2]
  out <- numerator/denominator
  if((!is.complex(q)) && !is.complex(a) && !is.complex(alpha)){ out <- Re(out)}
  return(out)
}

#' @export
`poch` <- function(a,q=a,n=Inf,maxit=1e6){ # n integer or infinite
  if(any(n != round(n))){return(poch_complex_n(a,q,n))}

  neg_found <- any(n<0)
  if(neg_found){
    wanted <- which(n<0)
    n <- abs(n)
    a[wanted] <- a[wanted]/q[wanted]^n[wanted]
  }

  jj <- cbind(Re(c(a)),Im(c(a)),Re(c(q)),Im(c(q)),c(n))
  jj <- pochhammer(a_real=jj[,1],a_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], n=jj[,5],maxit)
  
  z <- jj[,1] + 1i*jj[,2]
  if(neg_found){ z[wanted] <- 1/z[wanted] }

  if((!is.complex(q)) && !is.complex(a)){ z <- Re(z) }
  
  if(length(a) >= length(q)){
    attr <- attributes(a)
  } else {
    attr <- attributes(q)
  }

  if((length(a) == length(q)) && (length(a)==1) && (length(n)>1)){ attr <- attributes(n) }
  attributes(z) <- attr
  return(z)
}

#' @export
`pochprod` <- function(avec,q=1,n=Inf){ prod(sapply(avec,function(a){poch(a,q,n)}))}

#' @export
`pochR` <- function(a,q=a,n=Inf){
  if(any(n<0)){
    wanted <- which(n<0)
    n <- abs(n)
    q[wanted] <- a[wanted]/q[wanted]^n[wanted]
  }

  out <- 1
  outold <- 2
  subt <- a
  i <- 1
  while(any(outold != out) & (i<=n)){
    outold <- out
    out <- out*(1-subt)
    subt <- subt*q
    i <- i+1
  }
  
  if(any(n<0)){ out[wanted] <- 1/out[wanted] }
  return(out)
}

#' @export
rhs <- function(a,q,x){
  out <- 0
  outold <- 1
  n <- 0
  while(out != outold){
    outold <- out
    out <- out + poch(a,q,n)/poch(q,q,n) * x^n
    n <- n+1
  }
  return(out)
}
  

#' @export
disc <- function(a,x,q){

  LHS <- poch(a*x,q)/poch(x,q)
  RHS <- rhs(a,q,x)
  return(c(LHS=LHS,RHS=RHS,diff=RHS-LHS))
}

#' @export
q_factorial <- function(n,q=1){
    out <- poch(q,q,n)/(1-q)^n
    out[q==1] <- factorial(n[q==1])
    return(out)
}

#' @export
qfactorial_int <- function(n,q=1){
    stopifnot(n>=0)
    stopifnot(n==round(n))
    prod(unlist(lapply(sapply(seq_len(n),seq_len),\(x){sum(q^(x-1))})))
}

#' @export
q_choose <- function(n,k,q=1){q_factorial(n,q)/(q_factorial(n-k,q)*q_factorial(k,q))}

#' @export
q_gamma <- function(x,q=1){q_factorial(x-1,q)}

#' @export
q_binom <- function(n,k,q=1){ poch(q,q,n)/poch(q,q,k)/poch(q,q,n-k)}

#' @export
q_myexp <- function(z,q=1){
    jj <- cbind(c(Re(z)),c(Im(z)),c(Re(q)),c(Im(q)))
    jj <- qexp_C(z_real=jj[,1],z_imag=jj[,2], q_real=jj[,3],q_imag=jj[,4], maxit=1e6)
    out <- (jj[,1] + 1i*jj[,2])
    if((!is.complex(z)) && !is.complex(q)){z <- Re(z)}
    if(length(z) >= length(q)){attributes(out) <- attributes(z)} else {attributes(z) <- attributes(q)}
    return(out)
}

#' @export
q_exp <- function(x,q){1/poch(x,q)}

#' @export
q_Exp <- function(x,q){poch(-x,q)}

#' @export
q_sin <- function(x,q){(q_exp(1i*x,q) - q_exp(1i*x))/(2i)}

#' @export
q_cos <- function(x,q){(q_exp(1i*x,q) + q_exp(1i*x))/(2 )}

#' @export
q_Sin <- function(x,q){(q_Exp(1i*x,q) - q_Exp(1i*x))/(2i)}

#' @export
q_Cos <- function(x,q){(q_Exp(1i*x,q) + q_Exp(1i*x))/(2 )}

#' @export
q_beta <- function(a,b,q){q_gamma(a,q)*q_gamma(b,q)/qgamma(a+b,q)}

#' @export
q_n <- function(n,q=1){
    out <- (1-q^n)/(1-q)
    out[q==1] <- n
    return(out)
}

#' @export
q_hypergeo <- function(z,q, a,b, maxit = 2000){
  jj <- cbind(c(z),c(q))

  jj  <- BasicHypergeometric(z_real=Re(jj[,1]), z_imag=Im(jj[,1]), q_real=Re(jj[,2]), q_imag=Im(jj[,2]),
                                       a_real=Re(a), a_imag=Im(a), b_real=Re(b), b_imag=Im(b), maxit=maxit)
  return(jj)
}


