##WorkGroup 67; Members: Jonathan Hoover, Tianshu Liu, Zhiqing Wang
##GitRepo: 
##Practical 4: BFGS Optimization

##Goal:

#################################################################################################
##Start of Functions

##finite differentiation for gradient
fd_grad <- function(theta, f, fk, eps = 1e-7, ...){
  ##f: function to find the derivative for given parameters theta
  ##theta: values of parameters to evaluate f at
  ##fk: f evaluated at supplied values of theta
  ##eps: value close to zero to approximate a small change in "x"
  ##Finds derivative by: d/dtheta = (f(theta + eps) - f(theta))/eps
  ##for all theta_i in theta
  
  fd <- theta ##initialize fd with same size as theta
  
  for(i in 1:length(theta)){
    theta_i <- theta; theta_i[i] <- theta_i[i] + eps ##keeps all constant except for theta_i
    fi <- f(theta_i, ...) ##evaluate f at new value of theta with theta_i changed by eps
    fd[i] <- (fi - fk)/eps ##compute approximate derivative for theta_i
  }
  return(fd) ##return finite derivatives
}

fd_Hess <- function(theta, p, gradk, f, grad_attr, eps = 1e-7, ...){
  ##theta: values of parameters to evaluate f at
  ##p: number of parameters in theta
  ##gradk: gradient evaluated at values of theta (and f)
  ##f: objective function that gradient comes from; 
  ##   used to identify gradient
  ##grad_attr: Logical indicating if the gradient was supplied as a function
  ##eps: value close to zero to approximate a small change in "x"
  ##Finds derivative by: dgrad/dtheta = (grad(theta + eps) - grad(theta))/eps
  ##for all theta_i in theta
  
  Hfd <- matrix(0, nrow = p, ncol = p) ##initialize Hfd with zero pxp matrix
  
  ##get partial derivatives for all theta
  for(i in 1:p){
    theta_i <- theta; theta_i[i] <- theta_i[i] + eps ##keeps all constant except for theta_i
    
    if(grad_attr == T){
      fi <- f(theta_i, getg = T, ...) ##evaluate f at new theta with gradient
      gradi <- attr(fi, "gradient") ##pull from the gradient attr if included
    }else{
      fi <- f(theta_i, ...)##evaluate f at new theta without gradient
      gradi <- fd_grad(theta = theta_i, f = f, fk = fi, ...) ##calculate by finite differentiation if not
    }
    
    #Hfd <- append(Hfd,list((gradi - gradk)/eps)) ##compute approximate Hessian
    Hfd[i,] <- drop((gradi - gradk)/eps) ##compute approximate Hessian
  }
  
  return(Hfd) ##return finite derivatives
}


##bfgs
bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  ##

  p <- length(theta) ##number of parameters
  
  ####Initializing Values for loop and testing for finite values######
  
  I <- B0 <- diag(p) ##Set initial B (B0) to Identity matrix. I: used in computation of adjusted Hess
  theta_0 <- theta ##Current value of theta initialized with supplied theta
  f0 <- f(theta_0, ...) ##Evaluate function at initial theta
  grad_attr <- "gradient" %in% names(attributes(f0)) ##test if function evaluates with gradient
  
  ##Test if f evaluated at initial theta has a finite value
  ##Testing before gradient evaluation because it must be finite for gradient to work
  if(is.nan(f0) | is.infinite(f0)){
    stop("Evaluation of objective function at initial theta yielded a non-finite value")
  }
  
  ##set initial gradient based on whether gradient function is supplied with f
  if(grad_attr == T){
    grad0 <- attr(f0, "gradient") ##pull from the gradient attr if included
  }else{
    grad0 <- fd_grad(theta = theta_0, f = f, fk = f0, ...) ##calculate by finite differentiation if not
  }
  
  ##Test if initial Theta has a finite derivative
  if(sum(is.nan(grad0)) > 0 | sum(is.infinite(grad0)) > 0){
    stop(paste0("Evaluation of gradient at initial theta yielded a non-finite value\n", "f at initial theta: ", f0))
  }
   
  step <- -B0 %*% grad0 ##step defined by BFGS
  iter <- 0 ##count how many iterations we have done; start at 0
  c2 <- 0.9 ##Second wolfe condition c2
  step_iter <- 0 ##Count how many times we reduce step length (if too many issue a error)
  step_iter2 <- 0 ##Cont how many times we increase step length to meet wolf conditions
  
  ##Iterate to update theta, grad0, B, and f0 until the largest gradient value is smaller than
  ##the the function evaluated at the current theta * tolerance
  
  while(max(abs(grad0)) > (abs(f0) + fscale) * tol){
    theta_k <- theta_0 + step ##update theta to new theta_k
    fk <- f(theta_k, ...) ##evaluate function at new theta
    
    ##Test conditions for optimization
    ##1) New objective evaluation is finite
    ##2) New objective evaluation decreases relative to last evaluation
    ##3) Step Length satisfies the second wolf condition:
    ##   gradk^T * step >= c2 * grad0^T *step where c2 = 0.9
    
    ####Part 1: Finding Correct Step Length####
    
    ##reduce step length if non-finite or if step length doesn't decrease the function eval
    if(is.nan(fk) | is.infinite(fk) | (fk > f0)){
      
      step <- step * 0.5 ##reduce step length if it yielded a non-finite function eval
      step_iter <- step_iter + 1 ##count number of attempts at reducing step length
      
      if(step_iter > 1000){##1000 chosen arbitrarily to prevent an infinite (or very large loop) but to allow some testing
        stop(paste0("Attempts to reduce step size exceeded 1000 iterations before convergence at:\n", 
        "update iteration: ", iter, "\n",
        "f: ", fk, "; theta: ", paste(theta_k, collapse = ", "), "\n",
        "step length: ", paste(step, collapse = ", ")))
      }
      
      next ##start current iteration over again if it doesn't work
    }
    
    step_iter <- 0 ##reset step iteration for next update loop
    
    ##compute gradients for second wolf test
    if(grad_attr == T){
      gradk <- attr(fk, "gradient") ##pull from the gradient attr if included
    }else{
      gradk <- fd_grad(theta = theta_k, f = f, fk = fk, ...) ##calculate by finite differentiation if not
    }
    
    
    ##Test if new theta has a finite derivative (if not yield error)
    if(sum(is.nan(gradk)) > 0 | sum(is.infinite(gradk)) > 0){
      stop(paste0("Evaluation of gradient at theta_k yielded a non-finite value at:\n", 
                  "update iteration: ", iter, "\n",
                  "f: ", fk, "; theta: ", paste(theta_k, collapse = ", "), "\n"))
    }
    
    ##Once the value of fk < f0, test if it meets second wolfe condition
    ##If it does not, will need to increase slightly (since we reduced earlier)
    ##where c2 = 0.9
    if((t(gradk) %*% step) < (c2 * t(grad0) %*% step)){
      
      step <- step * 1.02 ## increase step slightly until the wolfe condition is met
      step_iter2 <- step_iter2 + 1 ##count number of attempts at increasing step_length to meet wolf condition
      
      if(step_iter2 > 1000){##1000 chosen arbitrarily to prevent an infinite (or very large loop) but to allow some testing
        stop(paste0("Attempts to meet second wolf condition exceeded 1000 iterations before convergence at:\n", 
                    "update iteration: ", iter, "\n",
                    "f: ", fk, "; theta: ", paste(theta_k, collapse = ", "), "\n",
                    "step length: ", paste(step, collapse = ", ")))
      }
      
      next ##repeat iteration with new step if it doesn't work until it does
    }
    
    step_iter2 <- 0 ##reset step_iteration (for the wolf conditions) for next update loop
    
    ####Part 2: Updating New Values####
    
    ##Calculate new B
    yk <- gradk - grad0; rho_k <- drop((t(step) %*% yk)^-1) ##rho_k is scalar
    
    ##Update for efficiency
    Bk <- (B0 - rho_k * step %*% (t(yk) %*% B0)) %*% (I - rho_k * yk %*% t(step)) + rho_k * step %*% t(step)
    # Bk <- (I - rho_k * step %*% t(yk))%*%B0%*%(I - rho_k * yk %*% t(step)) + rho_k * step %*% t(step)
    # Bk <- (B0 - rho_k * step %*% (t(yk) %*% B0)) ##intermediate step
    # Bk <- (Bk - (((Bk %*% yk) %*% t(step)) * rho_k)) - (rho_k * step %*% t(step)) ##final step
    
    ##Update B0 to B, thetak to theta0, fk to f0, etc once all conditions are met
    theta_0 <- theta_k; f0 <- fk; grad0 <- gradk;  B0 <- Bk;
    step <- -B0 %*% grad0 ##update step 0 to match new B
    iter <- iter + 1 ##update the iterations we have gone through (not counting step length adjustments)

    if(iter > maxit){
      stop("Surpassed max iterations before convergence")
    }
  }
  
  ##Compute approximate hessian at the last values of theta
  H <- fd_Hess(theta = drop(theta_0), p = p, gradk = grad0, f = f, grad_attr = grad_attr)
  
  ##Making symmetric for non-symmetric hessians derived with finite differentiation
  if(grad_attr == F){H <- 0.5 * (t(H) + H)}
                  
  return(list(f = f0[1], theta = drop(theta_0), g = drop(grad0), iter = iter, H = H))
    
}
##End of Functions
################################################################################################
theta <- c(-1,2)
rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by 'bfgs'
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb
f <- rb

bohachevsky <- function(theta, getg = F){
  x1 <- theta[1]; x2 <- theta[2]
  f <- x1^2 + 2*x2^2 - 0.3*cos(3*pi*x1 + 4*pi*x2) + 0.3
  if (getg) {
    attr(f, "gradient") <- g(theta) 
  }
  
  return(f)
}

beale <- function(theta){
  x <- theta[1]; y <- theta[2]
  f <- (1.5-x+x*y)^2+(2.25-x+(x*y)^2)^2 + (2.625-x+(x*y)^3)^2
  return(f)
}
# 
# ackley <- function(theta){
#   {\displaystyle f(x,y)=-20\exp \left[-0.2{\sqrt {0.5\left(x^{2}+y^{2}\right)}}\right]}
#   {\displaystyle -\exp \left[0.5\left(\cos 2\pi x+\cos 2\pi y\right)\right]+e+20}{\displaystyle -\exp \left[0.5\left(\cos 2\pi x+\cos 2\pi y\right)\right]+e+20}
# }



test <- bfgs(theta, rb, getg = T)
test2 <- bfgs(theta, rb)
optim(theta, rb, method = "BFGS", hessian = T)

bfgs(c(100,100), f = beale)
optim(c(4,1), fn = beale, hessian = T)
#bfgs(theta = c(0,0), f= bohachevsky, getg = T)
