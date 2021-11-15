rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by âbfgsâ
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

difference <- function(theta,f,...){
  
  len_theta <- length(theta)
  diff_results <- vector("numeric", len_theta)
  for (i in 1:len_theta){
    new_theta <- theta
    new_theta[-i] <- theta[-i]
    new_theta[i] <- theta[i]+ 1e-8
    diff_results[i] <- (f(new_theta,...) - f(theta,...))/1e-8}
  grad <- diff_results
  return (grad)} 

gradients <- function(theta,f,...){
  obj <- f(theta,...)
  if (is.null(attr(obj,"gradient"))){
    return (difference(theta,f,...))
  } else{
    return (attr(obj,"gradient"))
  }
}

wolfcon1<-function(theta_t,delta_t,grad,alpha_t, c1 = 0.1,f,...){
  f(theta_t+alpha_t*delta_t,...) - f(theta_t,...) - alpha_t*c1*t(grad)%*%delta_t
}

wolfcon2<-function(theta_t,delta_t,grad,alpha_t, c2 = 0.9,f,...){
  gradients(theta_t+alpha_t*delta_t,f,...)%*%delta_t - c2*grad%*%delta_t
}


alpha_determine<- function(theta_t,delta_t,grad,alpha_t, c1=0.1,c2=0.9,f,...){
  c1 <- 0.1
  c2 <- 0.9
  a <- 0
  b <- 1000000
  alpha <- 1
  j = 0
  while(TRUE){
    if (wolfcon1(theta_t,delta_t,grad,alpha_t = alpha, c1 = 0.1 , f,...) > 0){
      j <- j + 1
      b <- alpha
      alpha <- (a+ alpha)/2
      #print(alpha)
      next
    }
    if (wolfcon2(theta_t,delta_t,grad,alpha_t = alpha, c2 = 0.9 , f,...) < 0){
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
  obj <- f(theta,...)
  #  if (is.null(attr(obj,"gradient"))) {
  #    grad  <- difference(theta,f,...)
  #  } else{
  #      grad <- attr(obj,"gradient")}
  grad <- gradients(theta,f,...)
  t <- 0
  p <- length(theta)
  theta_t <- theta
  B_t <- diag(p)
  g <- grad
  f0 <- obj
  graphs <- matrix(0,maxit,p)
  while (TRUE){
    #cat("t",t)
    #cat("g",grad) #grad <- c(0.5946788, -2.4306700)
    #cat("theta",theta_t) #theta_t <- c(0.3119817, -0.5312700)
    #cat("B",B_t)  # B_t <- c(0.28503980, 0.07089414,0.07089414, 0.01763255); B_t<- matrix(B_t,nrow=2)
    if (t == maxit){break} else if (max(abs(g)) < (abs(f0)+fscale)*tol) {break}
    delta_t <- -B_t%*%grad
    alpha <- alpha_determine(theta_t,delta_t,grad,alpha_t, c1 = 0.1, c2 = 0.9,f,...)
    print(alpha)
    s_t <- alpha*delta_t
    y_t <- gradients(theta+ s_t,f,...) - gradients(theta,f,...)
    rho_t <- ((t(s_t)%*%y_t)[1,1])^(-1)
    sy_t <- s_t%*%t(y_t)
    
    t <- t+1
    #B_t <- B_t -rho_t*(B_t%*%t(sy_t)) - (rho_t)*((sy_t)%*%B_t) + (rho_t^2)*sy_t%*%B_t%*%t(sy_t) + rho_t*s_t%*%t(s_t)
    B_t <- (diag(p) - rho_t*sy_t)%*%B_t%*%t(diag(p) - rho_t*sy_t) + rho_t*s_t%*%t(s_t)
    theta_t <- theta_t + s_t
    grad <- gradients(theta_t,f,...)
    g <-grad
    obj <- f(theta_t,...)
    f0 <- obj
    graphs[t,] <- theta_t
  }
  return (graphs)
}

theta <- c(-2,1)
graphs <- bfgs(theta,rb,getg = TRUE, tol=1e-5,fscale=1,maxit=1000)  
graphs