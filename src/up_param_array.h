void up_vpar(const double rho, const double rho2,
             const arma::field<arma::uvec> uid,
             arma::field<arma::mat> & alpha,
             arma::field<arma::mat> & alpha_s,
             arma::mat & beta, arma::mat & beta_s);

void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta, 
                const int & k);