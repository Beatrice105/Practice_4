# Practice 4
## Aim: to write an R function, bfgs, implementing 
##      the BFGS quasi-Newton minimization method.
## Considerations:
## Outline high level design:

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
  if (f not supply gredient) {
    g <- ## compute the gradient by finite differencing
  }
  
  B <- diag(p)
  I <- diag(p)
  g <- f() ## 求法？？？
  theta.set <- matrix(0, 1, 2)
  theta.set[1, ] <- theta
  iter <- 1 ## initialize the number of iterations
  
  while (err0 > maxit) {
    
    ## the process of Newton 
    ## iter = 1
    theta <- theta.set[iter, ]
    f <- f(theta)
    g <- g(theta)
    
    
    s <- 
      theta.set[iter + 1, ] <- theta - s  ## theata(k+1) = theata(k) - B(k)g(k)
    y <- g(theta.set[iter + 1, ]) - g(theta)
    
    rho <- drop( solve(t(s) %*% y) ) ## simplify
    item <- (I - rho * s %*% t(y)) ## avoid repeat computation
    B <- item %*% B %*% item + rho * s %*% t(s)
    
    
    ## 迭代终止条件1：达到最大迭代次数
    if (iter == iterMax) {
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
rb <- function(theta,getg=FALSE,k=10) {
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















