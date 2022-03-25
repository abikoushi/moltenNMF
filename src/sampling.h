
arma::vec r_norm_vec(arma::vec mean, arma::vec sd);

void sample_mu(arma::mat & mu, arma::vec y, int N, arma::uvec xi, arma::uvec xp, arma::uvec varind, int K, int L, double lambda, double tau);

void up_ms(arma::mat & m, arma::mat & s2, arma::vec y, int N, arma::uvec xi, arma::uvec xp, arma::uvec varind, int K, int L, double lambda, double tau);
