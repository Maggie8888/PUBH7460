eliptical <- function(loglike,Y, Sig, N_mcmc,burn_in,N,flag){
  mcmcSamples <- matrix(rep(0,(burn_in+N_mcmc)*N),nrow=N_mcmc+burn_in,byrow=T)
  if(flag==TRUE){
    normSamples <- mvrnorm(n = burn_in+N_mcmc, mu = rep(0,N), Sigma = Sig , tol = 1e-6)
  }
  else{
    normSamples <- rcpp_rmvnorm(n = burn_in+N_mcmc, S = Sig, mu = rep(0,N))
  }
  unifSamples <- runif(n=burn_in+N_mcmc)
  theta <- runif(n=N_mcmc+burn_in,min=0,max=2*pi)
  theta_min <- theta-2*pi
  theta_max <- theta+2*pi
  for(i in 2:(N_mcmc+burn_in)){
    f <- mcmcSamples[i-1,]
    llh_thresh <- loglike(f,Y) + log(unifSamples[i])
    f_star <- f*cos(theta[i])+normSamples[i,]*sin(theta[i])
    while(loglike(f_star,Y) < llh_thresh)
    {
      if (theta[i] < 0) 
      {
        theta_min[i] <- theta[i]
      }
      else
      {
        theta_max[i] <- theta[i]
      } 
      theta[i] <- runif(n=1,min=theta_min[i],max=theta_max[i])  
      f_star <- f*cos(theta[i])+normSamples[i,]*sin(theta[i])
    }
    mcmcSamples[i,] <- f_star
  }
  return(mcmcSamples[(burn_in+1):(burn_in+N_mcmc),])
}
