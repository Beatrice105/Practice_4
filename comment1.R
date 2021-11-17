wolfcon1<-function(theta,delta,grad,alpha, c1 = 0.1,f,...){
  ### the function of wolf condition 1 ###
  f(theta+alpha*delta,...) - f(theta,...) - alpha*c1*t(grad)%*%delta
}

wolfcon2<-function(theta,delta,grad,alpha, c2 = 0.9,f,...){
  ### the function of wolf condition 2 ###
  t(gradients(theta+alpha*delta,f,...))%*%delta - c2*t(grad)%*%delta
}

alpha_determine<- function(theta,delta,grad, c1=0.1,c2=0.9,f,...){
  ### the general idea of this function is in paper at page5 "Dong Y. Step lengths in BFGS method for monotone gradients. Comput Math Appl. 2010; 60(3): 563- 571"
  ### our code also refer the code https://blog.csdn.net/qq_44401469/article/details/106392675
  ### the function for determining the reasonable value of alpha ###
  ### the basic idea of function to determine alpha is to make alpha both fulfill the wolfcondition1 and wolfcondition2
  ### if the theta is not fulfill the wolf condition1, we need to decrease alpha so that it could fulfill 
  ### if the theta is not fulfill the wolf condition2, we need to increase alpha so that it could fulfill
  
  ### The directly general idea is if theta disobey wolf1, then alpha = 0.5*alpha, if theta disobey wolf2, then alpha*2*alpha
  ### However, the question is that 
  ### "if the function is not fulfill the wolf condition1 at last iteration, then we half the alpha so that it could fulfill the wolfcon1
  ### but if the alpha now did not fulfill the wolf condition2, and if we still double the alpha now, we will go back original alpha"
  ### To avoid the above situations, we would record the last iteration alpha value.
  
  ### b is to record the alpha's value at the last term (the term which is most approach to now) where alpha disobey the wolf 1
  ### a is to record the alpha's value at the last term (the term which is most approach to now) where alpha disobey the wolf 2
  ### Now if the alpha's value disobey the wolf2, we will make alpha[t] = min(2*alpha[t-1], 0.5*(alpha[t-1] + b))
  ### the only 2*alpha[t-1] may make the alpha change from 1->0.5->1, the 0.5*(alpha[t-1] + b) could make alpha change from 1->0.5->0.75
  ### In other words, we could increase the alpha value at wolfcon2 but the value of alpha after increase should still smaller than the last term alpha value which disobey the wolfcon1
  ### Similarly, we could decrease the alpha value at wolfcon1 but the value of alpha after decrease should still larger than the last term the value disobey the worlf2
  ### The update alpha = 0.5*(a+alpha) would avoid 1->2->1, and make it become 1->2->1.5
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