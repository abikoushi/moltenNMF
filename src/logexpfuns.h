double xlogy(const double & x, const double & y);

double logsumexp(const arma::rowvec & x);

arma::rowvec softmax(const arma::rowvec & x);

arma::mat mat_digamma(const arma::mat & a);

arma::vec vec_digamma(arma::vec a);

void up_log_gamma(arma::mat & logv, const arma::vec & a, const double & logb, const int & l);

void up_log_gamma_uid(arma::mat & logv, const arma::vec & a, const double & logb,
                      const int & l, const arma::uvec uid);

double kld2(const arma::field<arma::mat> & alpha,
            const arma::mat & beta,
            const double & a,
            const double & b);
