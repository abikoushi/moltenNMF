void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta);

void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta, 
                const int & k);


void up_latentV_uid(arma::field<arma::mat> & V,
                    arma::field<arma::mat> & logV,
                    const arma::field<arma::mat> & alpha,
                    const arma::mat & beta, 
                    const int & k,
                    const arma::field<arma::uvec> & uid);

void up_vpar(const double rho,
             const double & rho2,
             arma::field<arma::mat> & alpha,
             const arma::field<arma::mat> & alpha_s,
             arma::mat & beta,
             const arma::mat & beta_s);

void up_vpar(const double rho,
             const double & rho2,
             const arma::field<arma::uvec> & uid,
             arma::field<arma::mat> & alpha,
             const arma::field<arma::mat> & alpha_s,
             arma::mat & beta,
             const arma::mat & beta_s);

void up_vpar(arma::field<arma::mat> & alpha,
             arma::mat & beta,
             arma::field<arma::mat> & V,
             arma::field<arma::mat> & logV,
             const arma::field<arma::mat> & alpha_s,
             const arma::mat & beta_s,
             const arma::field<arma::uvec> & uid,
             const double rho,
             const double & rho2);

