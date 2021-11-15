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
  theta.set <- matrix(0, 1, 2)
  theta.set[1, ] <- theta
  iter <- 1 ## initialize the number of iterations
  
  while (err0 > maxit) {
    
    ## the process of Newton 
    ## iter = 1
    theta <- theta.set[iter, ]
    f <- f(theta)
    g <- g(theta)
    
    ## Wolfe condition:
    alpha <- 1 ## initialize alpha = 1
    c1 <- 
    c2 <- 0.9
    delta <- - B %*% g(theta)
    f(theta + alpha * delta) <= f(theta) + c1 * t(g(theta)) %*% delta
    
    t(g(theta + alpha * delta)) >= c2 * t(g(theta)) %*% delta
   
    s <- 
      theta.set[iter + 1, ] <- theta - s  ## theata(k+1) = theata(k) - B(k)g(k)
    y <- g(theta.set[iter + 1, ]) - g(theta)
    
    rho <- drop( solve(t(s) %*% y) ) ## simplify
    item <- (I - rho * s %*% t(y)) ## avoid repeat computation
    B <- item %*% B %*% item + rho * s %*% t(s)
    
    
    ## 迭代终止条件1：达到最大迭代次数
    if (iter == maxit) {
      print("Reach the maximum number of BFGS iterations!")
      break
    }
    
    ## 迭代终止条件2：找到满足精度要求的解
    if (abs(theta(iter+1)-theta(iter)) < tol) {
      return()
      break
    }
    
    iter <-  iter + 1 
    
  }
  
  
  
  result <- c(f, theta, iter, g, H)
  return(result)
  
} 

theta <-  c(-1,2)
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


## theta如果是多维的，应该同时分别算？？？
theta ## 初始值
err0 <- inf
tol ## 最大容许误差
maxit ## 最大迭代次数















