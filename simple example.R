library(mvnfast)
p <- 2
y <- c(2, 0)
mu_0 <- c(-4, 5)
Sigma_0 <- diag(2)
chol_Sigma_0 <- chol(Sigma_0)
Sigma <- diag(2)
chol_Sigma <- chol(Sigma)
## Analytic Posterior
Sigma_post <- solve(solve(Sigma) + solve(Sigma_0))
mu_post <- Sigma_post %*% (solve(Sigma) %*% y + solve(Sigma_0) %*% mu_0)
## number of replicates for plotting
N <- 1000
x1 <- seq(-10, 10, length.out=200)
x2 <- seq(-10, 10, length.out=200)
z_prior <- matrix(dmvn(as.matrix(expand.grid(x1, x2)), mu_0, Sigma_0), length(x1), length(x2))
z_likelihood <- matrix(dmvn(as.matrix(expand.grid(x1, x2)), y, Sigma), length(x1), length(x2))
z_post <- matrix(dmvn(as.matrix(expand.grid(x1, x2)), mu_post, Sigma_post), length(x1), length(x2))
## Make Plot
contour(x1, x2, z_prior, col="blue")
contour(x1, x2, z_likelihood, col="red", add=TRUE)
contour(x1, x2, z_post, col="black", add=TRUE)
## log likelihood function
elliptical_slice <- function (mu, mu_mean, prior) { #, angle_range, varargin
        ## mu is the current value of the parameter
        ## prior is either a prior sample or a Cholesky of prior covariance
        # p <- length(mu)
        angle_range <- 0
        current_log_like <- dmvn(mu, y, Sigma, log=TRUE)
        # if (length(prior) == p) {
        ## user provided a prior sample
        nu <- prior
        # } else {
        #   ## user provided a Cholesky of prior covariance
        # nu <- t(prior) %*% rnorm(p)
        # }
        
        hh <- log(runif(1)) + current_log_like
        
        ## Setup a bracket and pick a first proposal
        ## Bracket whole ellipse with both edges at first proposed point
        phi <- runif(1) * 2 * pi
        phi_min <- phi - 2 * pi
        phi_max <- phi
        
        test <- TRUE
        ## Slice sampling loop
        while (test) {
                ## compute mu_proposal for proposed angle difference and check to see if it is on ther slice
                ## adjust mu so that it is mean 0
                mu_proposal <- (mu - mu_mean) * cos(phi) + nu * sin(phi)
                ## adjust the proposal so that it has the correct mean
                mu_proposal <- mu_proposal + mu_mean
                
                current_log_like <- dmvn(mu_proposal, y, Sigma, log=TRUE)
                if (current_log_like > hh) {
                        ## proposal is on the slice
                        mu <- mu_proposal
                        test <- FALSE
                } else if (phi > 0) {
                        phi_max <- phi
                } else if (phi < 0) {
                        phi_min <- phi
                } else {
                        warning("Bug detected - elliptical slice sampling shrunk to current position and still not acceptable")
                }
                ## Propose new angle difference
                phi <- runif(1) * (phi_max - phi_min) + phi_min
        }
        return(mu)
}

## Sample using MCMC  
n_mcmc <- 10000
mu <- matrix(0, n_mcmc, p)
mu[1, ] <- rmvn(1, rep(0, p), Sigma_0) + mean(mu)
for (i in 2:n_mcmc) {
        if (i %% 1000 == 0){
                cat(i, "\n")
        }
        mu_prior <- rmvn(1, rep(0, p), Sigma_0)
        mu[i, ] <- elliptical_slice(mu[i-1, ], mu_0, mu_prior)
}
## Make Plot
contour(x1, x2, z_prior, col="blue")
contour(x1, x2, z_likelihood, col="red", add=TRUE)
contour(x1, x2, z_post, col="black", add=TRUE)
points(mu, col=adjustcolor('black', alpha.f = 0.05))

Rcpp::sourceCpp("./ess.cpp")
## Sample using MCMC  
n_mcmc <- 10000
mu <- matrix(0, n_mcmc, p)
mu[1, ] <- rmvn(1, rep(0, p), Sigma_0)
for (i in 2:n_mcmc) {
        if (i %% 1000 == 0){
                cat(i, "\n")
        }
        mu_prior <- rmvn(1, rep(0, p), Sigma_0)
        mu[i, ] <- essVec(mu[i-1, ], mu_prior, mu_0, y, Sigma)
}
## Make Plot
contour(x1, x2, z_prior, col="blue")
contour(x1, x2, z_likelihood, col="red", add=TRUE)
contour(x1, x2, z_post, col="black", add=TRUE)
points(sweep(mu, 2, mu_0, "+"), col=adjustcolor('black', alpha.f = 0.05))