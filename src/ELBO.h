double lowerbound_logML_pois(const arma::mat & alpha,
                             const arma::mat & beta,
                             const arma::mat & lambda,
                             const arma::mat & loglambda,
                             const arma::vec & R,
                             const arma::vec & y,
                             const double & a,
                             const double & b);

double lowerbound_logML_pois_sp(const arma::mat & alpha,
                                const arma::mat & beta,
                                const arma::mat & lambda,
                                const arma::mat & loglambda,
                                const arma::vec & R,
                                const arma::vec & yv,
                                const arma::uvec & yi,
                                const double & a,
                                const double & b);
