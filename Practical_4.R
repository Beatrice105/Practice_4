# Practice 4
## Aim: to write an R function, bfgs, implementing 
##      the BFGS quasi-Newton minimization method.
## Considerations:
## Outline high level design:

grad <- function(theta, f, ..., eps = 1e-7) {
  f.value <- f(theta, ...)
  p <- length(theta)
  grad <- vector("numeric", p)
  if (is.null(attr(f.value,"gradient"))) { ## if f doesn't supply a gradient attribute, we compute the gradient by finite differencing
    for (i in 1:p) {
      theta.new <- theta
      theta.new[i] <- theta[i] + eps
      grad[i] <- (f(theta.new) - f(theta)) / eps
    }
    return(grad)
  } else { ## if f supplies a gradient attribute
    return(attr(f.value, "gradient"))
  }
}
  
Hessian <- function(theta, f, ..., eps = 1e-7) {
  f.value <- f(theta, ...)
  p <- length(theta)
  Hessian <- matrix(0, p, p)
  if (is.null(attr(f.value,"hessian"))) { ## if f doesn't supply a hessian attribute, we compute the gradient by finite differencing
    for (i in 1:p) {
      theta.new <- theta
      theta.new[i] <- theta[i] + eps
      Hessian[i,] <- (grad(theta.new, f) - grad(theta, f)) / eps
    }
    return(Hessian)
  } else { ## if f supplies a gradient attribute
    return(attr(f.value, "hessian"))
  }
}

wolfcon1<-function(theta, delta, g, alpha, c1 = 0.1, f){
  ### the function of wolf condition 1 ###
  f(theta + alpha * delta) <= f(theta) + alpha * c1 * t(g) %*% delta 
}

wolfcon2<-function(theta,delta,g,alpha, c2 = 0.9,f){
  ### the function of wolf condition 2 ###
  t(grad(theta+alpha*delta,f)) %*% delta >= c2 * t(g) %*% delta
}

alpha_determine <- function(theta, delta, grad, c1=0.1, c2=0.9, f){
  ### the function for determining the reasonable value of alpha ###
  a <- 0
  b <- 1000000
  alpha <- 1
  j = 0
  while(TRUE){
    if (wolfcon1(theta,delta,grad,alpha, c1 , f)){
      ### For the value the alpha get now, if it disobey the wolfcon1, then execute the following ###
      j <- j + 1
      b <- alpha
      alpha <- (a+ alpha)/2
      #print(alpha)
      next
    }
    if (wolfcon2(theta,delta,grad,alpha, c2, f)){
      ### For the value the alpha get now, if it disobey the wolfcon2, then execute the following ###
      a <- alpha
      alpha <- min(c(2*alpha,(alpha+b)/2))
      #print(alpha)
      next
    }
    break
  }
  return (alpha)
}

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {
  
  ## input: 
  ##    theta: a vector of initial values for the optimization parameters.
  ##    f: the objective function to minimize. 
  ##    tol: the convergence tolerance.
  ##    fscale: a rough estimate of the magnitude of f at the optimum - used in convergence testing.
  ##    maxit: the maximum number of BFGS iterations to try before giving up.
  ## output: f, theta, iter, g, H
  ##    f: the scalar value of the objective function at the minimum. 
  ##    theta: the vector of values of the parameters at the minimum.
  ##    iter: the number of iterations taken to reach the minimum.
  ##    g: the gradient vector at the minimum (so the user can judge closeness to numerical zero).
  ##    H: the approximate Hessian matrix (obtained by ﬁnite diﬀerencing) at the minimum.
  
  p <- length(theta) ## number of parameters
  B <- diag(p) ## initialize 
  I <- diag(p)
  theta.set <- matrix(0, maxit, 2)
  theta.set[1, ] <- theta
  iter <- 0 ## initialize the number of iterations
  f0 <- f(theta)
  g <- grad(theta, f)
  
  while (max(abs(g)) >= (abs(f0)+fscale) * tol) {
    
    ## the process of Newton 
    iter <-  iter + 1
    theta <- theta.set[iter, ]
    g <- grad(theta, f)
    delta <- - B %*% g
    
    alpha <- alpha_determine(theta, delta, g, c1 = 0.1, c2 = 0.9, f)
    s <- alpha * delta
    theta.set[iter + 1, ] <- theta - s  ## theta(k+1) = theta(k) - s(k)
    y <- grad(theta.set[iter + 1, ], f) - g
    
    rho <- 1 / drop(t(s) %*% y)
    item <- I - rho * s %*% t(y) ## avoid repeat computation
    B <- item %*% B %*% item + rho * s %*% t(s)
    if (iter == maxit) break
    
  }
 
  g <- grad(theta, f)
  H <- Hessian(theta, f)
  
  result <- list(f = f0, theta = theta , iter = iter, g = g, H = H)
  return(result)
  
} 

theta <- c(-1,2)
f <- rb
rb <- function(theta, getg=FALSE, k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2), -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb
optim(c(-1,2), rb, hessian = T, method="BFGS")
bfgs(theta = c(-1,2), f = rb)

boha1 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- x1^2
  term2 <- 2*x2^2
  term3 <- -0.3 * cos(3*pi*x1)
  term4 <- -0.4 * cos(4*pi*x2)
  
  y <- term1 + term2 + term3 + term4 + 0.7
  return(y)
}

boha1(theta)
bfgs(theta, boha1)


