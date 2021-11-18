##WorkGroup 67; Members: Zhiqing Wang, Tianshu Liu, Jonathan Hoover, 
##GitRepo: https://github.com/Beatrice105/Practice_4
##Practical 4: BFGS Optimization

##Goal:
##The Overall goal of this project is to write a function to implement
##the quasi-newton BFGS optimization method
##The method attempts to minimize an objective function f given
##Starting parameters theta. The function uses successive approximations
##Of the inverse hessian to update a step length for theta until
##the gradient of the function evaluated at the new theta is approximately 0.
##The update of B (Inverse hessian) is defined such that the step length
##Meets the wolfe conditions: 
## 1) f(new theta) <= (f(old theta) + c1 * grad(old theta)^T * steplength)
## 2) (grad(new theta)^T * alpha * steplength)  >= (c2 * grad(old theta)^T * alpha*steplength)

##The following functions break up these necessary tasks
##But the main function is the last: bfgs

#################################################################################################
##Start of Functions


difference <- function(theta,f,eps = 1e-8,...){
  ##Goal: Perform finite differentiation on an objective function f
  ##      to find the derivatives (gradient) of f at theta. 
  ##      The idea is to perform (f(theta + eps) - f(theta))/eps
  ##      To approximate differentiation with a small value eps
  
  ##Inputs:
  ##theta: values of theta where f is to be differentiated
  ##f: Objective function to differentiate
  ##eps: small value used to approximate 0 when differentiating
  
  
  p <- length(theta) #number of parameters in theta
  diff_results <- vector("numeric", p) #create the empty parameter vector
  for (i in 1:p){##finite differentiation for each parameter to get gradient
    new_theta <- theta # copy the original theta to create new vector of theta
    new_theta[i] <- theta[i] + eps # replace the ith parameter with the theta[i] + eps
    diff_results[i] <- (f(new_theta,...) - f(theta,...))/eps ##f(x_new) - f(x_previous)/ eps
    } 
  # output the gradient
  return (diff_results)} 


gradients <- function(theta,f,...){
  ##Goal:
  ##function to calculate gradients whether gradient functions are supplied or not
  
  ##Inputs:
  ##theta: parameter values to evaluate gradient at
  ##f: objective function to evaluate at theta (for finite differentiation step)
  
  ##determine way to calculate gradient
  obj <- f(theta,...)
  if (is.null(attr(obj,"gradient"))){ # if there does not exist the function to give the gradient function 
    return (difference(theta = theta,f = f,...)) # use difference function to calculate the derivative
  } else{ # if there exists the function to give the gradient
    return (attr(obj,"gradient")) # output the gradient value
  }
}

fd_Hess <- function(theta, gradk, f, eps = 1e-7,...){
  ##Goal: function to evaluate Hessian of a given function f at values theta
  ##      via finite differentiation of gradients
  ##      Finds derivative by: dgrad/dtheta = (grad(theta + eps) - grad(theta))/eps

  
  ##theta: values of parameters to evaluate f at
  ##gradk: gradient evaluated at values of theta (and f)
  ##f: objective function that gradient comes from; 
  ##   used to identify gradient
  ##eps: value close to zero to approximate a small change in "x"
  ##...: Additional arguments of provided function
  
  
  p <- length(theta) ##number of parameters in theta
  Hfd <- matrix(0, nrow = p, ncol = p) ##initialize Hfd with zero pxp matrix
  
  ##get partial derivatives for all theta
  for(i in 1:p){
    theta_i <- theta; theta_i[i] <- theta_i[i] + eps ##keeps all constant except for theta_i
    gradi <- gradients(theta = theta_i, f=f, ...)
    
    Hfd[i,] <- drop((gradi - gradk)/eps) ##compute approximate Hessian; need to use drop to get correct dims
  }
  
  return(Hfd) ##return finite derivatives
}


wolfcon1<-function(f0,fk,alpha,delta,grad,c1 = 0.1,...){
  ##Goal: Test wolfe condition 1 (where 0 < c1 < c2 < 1)
  
  ##It evaluates the difference between left hand of wolfe condition 1
  ##And right hand of wolfe condition 1
  ##ie f(new theta) - (f(old theta) + c1 * grad(old theta)^T * steplength)
  ##ie if the result is <= 0  the wolfe condition is satisfied
  ##if the result is > 0 the wolfe condition is not satisfied
  
  ##Inputs:
  ##f0: objective function evaluated at previous theta
  ##fk: objective function evaluated at updated theta
  ##alpha: scaling factor for step length determined in another function
  ##delta: initial step length before reduction under wolfe condition 1
  ##grad: gradient evaluated at previous theta
  ##c1: Initial wolfe 1 condition c1 where  0 < c1 < c2 < 1; default 0.1

  
  #fk <- f(theta+alpha*delta,...) ##f at new theta
  return(fk - f0 - alpha*c1*t(grad)%*%delta)
}

wolfcon2<-function(alpha,delta,grad,gradk,c2 = 0.9,...){
  ##Goal: Test wolfe condition 2 (where 0 < c1 < c2 < 1)
  
  ##It evaluates the difference between left hand side
  ##and right hand side of wolfe condition 2 given a step length scaling by factor alpha
  ##ie (grad(new theta)^T * alpha * steplength) - (c2 * grad(old theta)^T * alpha*steplength)
  ##If the result is >= 0, wolfe condition 2 is satisfied
  ##if the result is < 0, wolfe condition 2 is not satisfied
  
  
  ##Inputs
  ##alpha: scaling factor of step length to ensure wolfe conditions are met
  ##delta: initial step length
  ##grad: gradient evaluated at un-updated theta values
  ##gradk: gradient evaluated at updated theta
  ##c2: Initial wolfe 2 condition c2 where  0 < c1 < c2 < 1; default 0.9
  
  #gradk <- gradients(theta+alpha*delta,f=f,...) ##gradient at updated theta
  return(t(gradk)%*%(alpha*delta) - c2*t(grad)%*%(alpha*delta)) ##test wolfe
}



alpha_determine<- function(theta,delta,grad, f, f0,iter, c1=0.1,c2=0.9,...){
  ##Goal:
  ##Evaluate a reasonable scaling factor for step length, called alpha here,
  ##To enable satisfaction of wolfe conditions
  
  ## The general idea of this function is in paper at page5 "Dong Y. 
  ## Step lengths in BFGS method for monotone gradients. Comput Math Appl. 2010; 60(3): 563- 571"
  ## our code also refer the code https://blog.csdn.net/qq_44401469/article/details/106392675
  
  ##Inputs:
  ##theta: previous value of theta
  ##delta: unscaled step length from BFGS algorithm
  ##grad: gradient evaluated at previous theta
  ##f: objection function to evaluate at new theta values (after consider alpha modification)
  ##f0: objective function evaluated at previous theta
  ##iter: iteration in bfgs loop that we are testing
  ##c1: Initial wolfe 1 condition c1 where  0 < c1 < c2 < 1; default 0.1
  ##c2: Initial wolfe 2 condition c2 where 0 < c1 < c2 < 1; default 0.9

  ### Methodology:
  ### the goal of the function is to determine a scaling factor, alpha, that lets step length fulfill the wolf conditions 1 and 2
  ### if the theta does not fulfill wolf condition1, we need to decrease alpha so that it can fulfill 
  ### if the theta does not fulfill wolf condition2, we need to increase alpha so that it can fulfill
  ### The remaining comments give a detailed explanation of the methodology; feel free to skip if this is sufficient.
  
  ### The idea is that if theta disobeys wolf1, then we initially reduce the step length until it fulfills wolf1.
  ### However, too much of a decrease may prevent it from fulfilling wolf2. So if the reduction does not meet wolf2,
  ### We must slightly increase alpha to meet wolf2. However, we still need to ensure that the subsequent increase does not surpass 
  ### the last iteration of alpha that failed wolf1. In other words, the increase still needs to fulfill wolf1.

  ### To ensure this, we record both the iteration of alpha that passed wolf1 and the last iteration of alpha that failed wolf1
  ### b is to record alpha's value at the last iteration where alpha did not pass wolf1
  ### a is to record alpha's value at the last iteration where alpha passed wolf1 but not wolf2
  
  ### If alpha disobeys wolf2, we update alpha to be alpha[t] = min(2*alpha[t-1], 0.5*(alpha[t-1] + b))
  ### In other words we take the smallest increase from either 2*(alpha that passes wolf1) or 
  ### mean(alpha that passed wolf 1, alpha that failed wolf 1)
  
  ### In essence, this strategy increases alpha to satisfy wolf2 but not so much that it should not satisfy wolf1
  ### We also need to ensure that any subsequent decrease during wolf1 does not decrease so much that it again does not satisfy wolfe 2
  ### This is why the first step updates alpha to be 0.5*(a + alpha) = mean(alpha that passes wolf 1, increased alpha)
  
  ### Again, these steps ensure that alpha gets more and more narrowed with each iteration
  ### We will refer to this methodology section throughout the code
  
  
  ### the function for determining the reasonable value of alpha ###
  alpha <- 1 ##initial scaling factor for step (ie no reduction)
  b <- Inf ##previous value of alpha that did NOT pass condition 1; Start large to ensure first step takes alpha = (2*alpha)
  a <- 0 ##value to scale alpha by (either 0 or prev alpha that DID pass cond 1)
  step_iter1 <- 0 ##initialize how many iterations have passed to reduce alpha (condition 1)
  step_iter2 <- 0 ##initialize how many iterations have passed to increase alpha (condition 2)
  
  while(TRUE){
    theta_k <- theta + alpha*delta ##theta updated with new alpha
    fk <- f(theta_k,...) ##f(theta) updated with new theta and new alpha
    
    ##Test to make sure new objective before wolfe 1 update is finite
    if(is.nan(fk) | is.infinite(fk)){
      alpha = alpha * 0.5 ##reduce alpha by half to make finite
      next ##start iteration again if evaluation of fk is non-finite
    }
    
    
    ##Wolf condition 1
    if (wolfcon1(f0=f0,fk = fk,alpha,delta,grad, c1 ,...) > 0){
      
      ##See methodology for a full descriptin of these steps
      
      step_iter1 <- step_iter1 + 1 ##iteration to count how many attempts at reduction
      b <- alpha ##store most recent alpha; used to update alpha for wolfe condition 2
      alpha <- (a+ alpha)/2 ##scale down using value for "a"; See notes above
      
      ##if too many attempts at reduction occur without a reduction  issue an error
      if(step_iter1 > 5000){
        stop(paste0("Attempts to reduce step size exceeded 5000 iterations before convergence at:\n",
                    "update iteration: ", iter, "\n",
                    "f: ", fk, "; theta: ", paste(theta_k, collapse = ", "), "\n",
                    "gradient: ", paste(gradients(theta + delta*alpha,f = f,...), collapse = ", "), "\n",
                    "step length: ", paste(delta*alpha, collapse = ", ")))
      }
      
      next ##attempt wolfe condition 1 again until met
    }
    
    step_iter1 <- 0 ##reset step iteration count for the first condition test
    gradk <- gradients(theta_k,f=f,...) ##gradients updated at new theta and alpha after wolf 1
    
    
    ##Test to make sure gradient is finite; if it is not, it will cause problems for wolfe condition 2
    if(sum(is.nan(gradk)) > 0 | sum(is.infinite(gradk)) > 0){
      alpha = alpha * 0.5 ##reduce alpha by half to make finite
      warning("Evaluation of gradient yielded a non-finite value")
      next ##start iteration again if evaluation of gradient is non-finite
    }
    
    
    ##Test Wolfe Condition 2
    if (wolfcon2(alpha,delta,grad,gradk = gradk,c2,...) < 0){
      
      ##See methodology above for detailed explanation of this algorithm
      step_iter2 <- step_iter2 +1
      a <- alpha
      alpha <- min(c(2*alpha,(alpha+b)/2))
      
      ##if too many attempts at increase occur without a proper increase issue an error
      if(step_iter2 > 5000){
        stop(paste0("Attempts to meet second wolf condition exceeded 5000 iterations before convergence at:\n",
                    "update iteration: ", iter, "\n",
                    "f: ", fk, "; theta: ", paste(theta_k, collapse = ", "), "\n",
                    "gradient: ", paste(gradk, collapse = ", "), "\n",
                    "step length: ", paste(delta*alpha, collapse = ", ")))
      }
      
      next##attempt wolfe condition 1 and 2 again until met
    }
    
    break ##once both are met exist the while looop
  }
  return (alpha) ##return the step length scaling factor that meets both wolfe conditions
}



bfgs <-function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  ##Goal:
  ##Write a function to implement the quasi-newton BFGS optimization method
  ##The method attempts to minimize an objective function f given
  ##Starting parameters theta. The function uses successive approximations
  ##Of the inverse hessian to update a step length for theta until
  ##the gradient of the function evaluated at the new theta is approximately 0.
  ##The update of B (Inverse hessian) is defined such that the step length
  ##Meets the wolfe conditions: 
  ## 1) f(new theta) <= (f(old theta) + c1 * grad(old theta)^T * steplength)
  ## 2) (grad(new theta)^T * alpha * steplength)  >= (c2 * grad(old theta)^T * alpha*steplength)
  
  ##Thus the first step is to find a scaling factor, alpha, such that the step length
  ##meets these conditions. And the second step to to update B
  ##We have defined 3 separate functions (wolfcon1, wolfcon2, and alpha_determine) above
  ##to solve the step length problem. Please refer to them for a detailed methodology of the process
  
  ##Inputs: 
  ##theta: initial value of parameters to start optimization process
  ##f: objective function to minimize (If the output of f provides a "gradient" attribute,
  ##   these values will be used as gradient; otherwise, the function will use finite differentiation to solve)
  ##...: additional parameters for f(...)
  ##tol: error tolerated between gradient and 0 (how small gradient must be to be considered optimized)
  ##fscale: scaling factor to account for objective functions that are optimized at (0,0)
  ##maxit: maximum number of iterations allowed to try to reach convergence before stopping the function
  
  ##Output:
  ##f: value of objective evaluated at the optimized theta values
  ##theta: value of parameters that optimizes f
  ##g: gradient of the objective (if f yields an object with a "gradient" attribute, it will solve via this;
  ##   otherwise it will solve via finite differentiation)
  ##iter: total number of iterations required to update theta to an optimized value
  ##H: Approximate Hessian at the optimized point solved via finite 
  ##   differentiation of gradients at optimized theta
  
  
  ####Initializing Values for update loop and testing for finite values######
  p <- length(theta) ##number of parameters in theta
  I <- B0 <- diag(p) ##Set initial B (B0) to Identity matrix. I: used in computation of adjusted Hess
  theta_0 <- theta ##Current value of theta initialized with supplied theta
  f0 <- f(theta_0, ...) ##Current evaluation of theta initialized with supplied theta
  
  
  ##Test if f evaluated at initial theta has a finite value
  ##Testing before gradient evaluation because it must be finite for gradient to work
  if(is.nan(f0) | is.infinite(f0)){
    stop("Evaluation of objective function at initial theta yielded a non-finite value")
  }
  
  ##set initial gradient for theta_t before updating (ie iteration 0)
  grad_0 <- gradients(theta = theta_0,f = f,...)
  
  ##Test if initial Theta has a finite derivative
  if(sum(is.nan(grad_0)) > 0 | sum(is.infinite(grad_0)) > 0){
    stop(paste0("Evaluation of gradient at initial theta yielded a non-finite value\n", "f at initial theta: ", f0))
  }
  
  iter <- 0 ##count how many iterations we have done; start at 0

  
  ##Iterate to update theta, grad0, B, and f0 until the largest gradient value is smaller than
  ##the the function evaluated at the current theta * tolerance
  
  while(max(abs(grad_0)) > (abs(f0) + fscale) * tol){
    
    if (iter > maxit){break} ##stop evaluating if convergence has not occurred before max # iterations
    
    delta <- -B0 %*% grad_0 ##delta at step iter = -B[iter]*gradient; ie prev -B * prev grad
    
    #determine the suitable alpha[iter] to scale the step to meet wolfe conditions
    ##see alpha_determine documentation for methodology
    alpha <- alpha_determine(theta_0,delta,grad_0, f = f, f0 = f0, 
                               c1 = 0.1, c2 = 0.9 , iter = iter,...) 
    
    ##calculate updated theta, step, and gradient with new step = s_t = alpha*delta
    ##where alpha scales delta to meet wolfe conditions
    
    s_t <- alpha*delta # scaled step to meet wolf conditions
    theta_t <- theta_0 + s_t ##updated theta at current iteration
    ft <- f(theta_t, ...) ##updated objective evaluated at new theta
    grad_t <- gradients(theta_t,f=f,...) ##updated gradient at current iteration
    
    ##Calculate values needed to update B
    y_t <- grad_t - grad_0 #difference between current grad_t and previous grad_0
    rho_t <- ((t(s_t)%*%y_t)[1,1])^(-1) #rho defined as (s_t^T*y_t)^-1; scalar

    ##Update B
    ##Normal formula for updating B is:
    #Bt <- (I - rho_t * s_t %*% t(y_t))%*%B0%*%(I - rho_t * y_t %*% t(s_t)) + rho_t * s_t %*% t(s_t)
    
    ##A more efficient formula is:
    update_coefficient <- (I - rho_t*(s_t%*%t(y_t))) ##multiply in this order for efficiency
    Bt <- (update_coefficient) %*% B0 %*% t(update_coefficient) + rho_t*s_t%*%t(s_t)
    
    ##Update theta_t to theta_0, ft to f0, and Bt to B0 etc once all conditions are met
    theta_0 <- theta_t; f0 <- ft; grad_0 <- grad_t; B0 <- Bt
    iter <- iter + 1 ##update the iterations we have gone through (not counting step length adjustments)

  }
  
  ##Calculate hessian by finite differentiation
  H <- fd_Hess(theta = theta_t, gradk = grad_t, f = f, eps = 1e-7, ...)
  
  ##Making symmetric for non-symmetric hessians derived with finite differentiation
  if(is.null(attr(f0, "gradient"))){H <- 0.5 * (t(H) + H)}
  
  ##return list of all necessary outputs
  return (list(f = f0[1], theta = drop(theta_t), iter = iter, g = drop(grad_t), H = H))
}


##End of Functions
#################################################################################################


nll3 <- function(theta,time,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (log(alpha),beta)
  alpha <- exp(theta[1]) ## so theta[1] unconstrained, but alpha > 0
  beta <- theta[2]
  mu <- alpha * exp(beta * time) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood
} #

t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases

bfgs(theta = c(log(10), 0.1), f = nll3, y = y, time = t80)
optim(c(log(10), 0.2), fn = nll3, y = y, time = t80, method = "BFGS", hessian = T)

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
#f <- rb

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


boha1 <- function(xx)
{
  ##########################################################################
  #
  # BOHACHEVSKY FUNCTION 1
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- x1^2
  term2 <- 2*x2^2
  term3 <- -0.3 * cos(3*pi*x1)
  term4 <- -0.4 * cos(4*pi*x2)
  
  y <- term1 + term2 + term3 + term4 + 0.7
  return(y)
}

ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  ##########################################################################
  #
  # ACKLEY FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a = constant (optional), with default value 20
  # b = constant (optional), with default value 0.2
  # c = constant (optional), with default value 2*pi
  #
  ##########################################################################
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}

camel6 <- function(xx)
{
  ##########################################################################
  #
  # SIX-HUMP CAMEL FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- (4-2.1*x1^2+(x1^4)/3) * x1^2
  term2 <- x1*x2
  term3 <- (-4+4*x2^2) * x2^2
  
  y <- term1 + term2 + term3
  return(y)
}



nll3 <- function(theta,time,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (log(alpha),beta)
  alpha <- exp(theta[1]) ## so theta[1] unconstrained, but alpha > 0
  beta <- theta[2]
  mu <- alpha * exp(beta * time) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood
} #

t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases



optim(c(log(10), 0.1), fn = nll3, y = y, time = t80, method = "BFGS", hessian = T)
bfgs(theta = c(log(10),0.1), f = nll3, y = y, time = t80)


bfgs(theta, rb, getg = T)
bfgs(theta, rb)
optim(theta, rb, method = "BFGS", hessian = T)

bfgs(c(4,1), f = beale)
optim(c(4,1), fn = beale, hessian = T)


bfgs(c(5,4), f = boha1)
optim(c(5,4), fn = boha1, method = "BFGS", hessian = T)

#bfgs(c(3,2), f = ackley)
optim(c(3,2), fn = ackley, method = "BFGS")

bfgs(c(10,10), f = camel6)
optim(c(10,10), f = camel6, method = "BFGS", hessian = T)

