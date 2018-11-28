// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Elliptical Slice Sampling

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double dMVNorm(const arma::vec & y, const arma::vec & mu, 
               const arma::mat & Sig, const bool logd = true){ 
        arma::mat rooti = trans(inv(trimatu(chol(Sig))));
        double rootisum = sum(log(rooti.diag()));      
        double constants = -(y.n_elem / 2.0) * log2pi;
        arma::vec z = rooti * (y - mu) ;      
        double out = constants - 0.5 * sum(z % z) + rootisum;       
        if(logd){
                return(out);
        } else {
                return(exp(out));
        }
}


// MCMC loop

// [[Rcpp::export]]
arma::vec essVec (const arma::vec& x_current, const arma::vec& x_prior, 
                  const arma::vec& mu_x, const arma::vec& y,
                  const arma::mat& Sigma) {
        // x_current is the current value of the parameter
        // x_prior is a prior sample
        // mu_x is the prior mean
        double current_log_like = dMVNorm(x_current + mu_x, y, Sigma, true);
        double hh = log(R::runif(0, 1)) + current_log_like;
        
        // Setup a bracket and pick a first proposal
        // Bracket whole ellipse with both edges at first proposed point
        double phi = R::runif(0, 1) * 2.0 * arma::datum::pi;
        double phi_min = phi - 2.0 * arma::datum::pi;
        double phi_max = phi;
        
        arma::vec x_ess = x_current;
        bool test = true;
        // Slice sampling loop
        while (test) {
                // compute mu_proposal for proposed angle difference and check to see if it is on ther slice
                arma::vec x_proposal = x_current * cos(phi) + x_prior * sin(phi);
                current_log_like = dMVNorm(x_proposal + mu_x, y, Sigma, true);
                if (current_log_like > hh) {
                        // proposal is on the slice
                        x_ess = x_proposal;
                        test = false;
                } else if (phi > 0) {
                        phi_max = phi;
                } else if (phi < 0) {
                        phi_min = phi;
                } else {
                        Rcout << "Bug detected - elliptical slice sampling shrunk to current position and still not acceptable";
                }
                // Propose new angle difference
                phi = R::runif(0, 1) * (phi_max - phi_min) + phi_min;
        }
        return(x_ess);
}
