rb <- function(theta,getg=FALSE,k=10) {
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

difference <- function(theta,f,...){
  ### the function for differentiation method to find derivative ###
  len_theta <- length(theta) # how many uncertain parameter theta
  diff_results <- vector("numeric", len_theta) # create the empty parameter vector
  for (i in 1:len_theta){
    new_theta <- theta # copy the original theta to create new vector of theta
    new_theta[i] <- theta[i]+ 1e-8 # replace the i th position of theta by the theta[i] +1e-8
    diff_results[i] <- (f(new_theta,...) - f(theta,...))/1e-8} # f(x_new) - f(x)/ 1e-8
  grad <- diff_results # output the gradient
  return (grad)} 


gradients <- function(theta,f,...){
  ### the function for getting the gradient of f
  obj <- f(theta,...) 
  if (is.null(attr(obj,"gradient"))){ # if there does not exist the function to give the gradient function 
    return (difference(theta,f,...)) # use difference function to calculate the derivative
  } else{ # if there exists the function to give the gradient
    return (attr(obj,"gradient")) # output the gradient value
  }
}

wolfcon1<-function(theta,delta,grad,alpha, c1 = 0.1,f,...){
  ### the function of wolf condition 1 ###
  f(theta+alpha*delta,...) - f(theta,...) - alpha*c1*t(grad)%*%delta
}

wolfcon2<-function(theta,delta,grad,alpha, c2 = 0.9,f,...){
  ### the function of wolf condition 2 ###
  t(gradients(theta+alpha*delta,f,...))%*%delta - c2*t(grad)%*%delta
}

alpha_determine<- function(theta,delta,grad, c1=0.1,c2=0.9,f,...){
  ### the function for determining the reasonable value of alpha ###
  a <- 0
  b <- 1000000
  alpha <- 1
  j = 0
  while(TRUE){
    if (wolfcon1(theta,delta,grad,alpha, c1 , f,...) > 0){
      ### For the value the alpha get now, if it disobey the wolfcon1, then execute the following ###
      j <- j + 1
      b <- alpha
      alpha <- (a+ alpha)/2
      #print(alpha)
      next
    }
    if (wolfcon2(theta,delta,grad,alpha, c2, f,...) < 0){
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

bfgs <-function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  t <- 0 # initial step t = 0
  obj <- f(theta,...) # get the objective value with theta
  grad_t <- gradients(theta,f,...) # the gradient value at step t = 0
  p <- length(theta) # length of theta 
  theta_t <- theta # the theta value at step t = 0
  B_t <- diag(p) # the initial inverse hensian matrix approximation B_t matrix at t= 0
  g <- grad_t 
  f0 <- obj 
  graphs <- matrix(0,maxit,p)
  while (TRUE){
    if (t == maxit){break} else if (max(abs(g)) < (abs(f0)+fscale)*tol) {break}
    delta_t <- -B_t%*%grad_t # delta at step t = - B[t]*gradient[t]
    alpha_t <- alpha_determine(theta_t,delta_t,grad_t, c1 = 0.1, c2 = 0.9,f,...) # determine the suitable alpha[t]
    s_t <- alpha_t*delta_t # s[t] = alpha[t]*delta[t]
    y_t <- gradients(theta_t+ s_t,f,...) - gradients(theta_t,f,...) # y[t] = gradient(theta[t]+s[t]) - gradient(theta[t])
    rho_t <- ((t(s_t)%*%y_t)[1,1])^(-1) # rho[t] = s[t]'*y[t]
    sy_t <- s_t%*%t(y_t)  #sy[t] = s[t]*y[t]'
    
    t <- t+1 # update step t= t+1
    #B_t <- B_t -rho_t*(B_t%*%t(sy_t)) - (rho_t)*((sy_t)%*%B_t) + (rho_t^2)*sy_t%*%B_t%*%t(sy_t) + rho_t*s_t%*%t(s_t)
    B_t <- (diag(p) - rho_t*sy_t)%*%B_t%*%t(diag(p) - rho_t*sy_t) + rho_t*s_t%*%t(s_t)
    theta_t <- theta_t + s_t # update the theta[t] = theta[t] + s[t]
    grad_t <- gradients(theta_t,f,...) # update gradient[t] = gradients(theta[t])
    g <-grad_t
    obj <- f(theta_t,...)
    f0 <- obj
    graphs[t,] <- theta_t
    print(alpha_t)
  }
  return (graphs)
}

theta <- c(-1,2)
graphs <- bfgs(theta,rb,tol=1e-5,fscale=1,maxit=100)  
graphs
