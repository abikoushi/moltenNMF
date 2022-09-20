arma::mat mysum_t(int n, arma::uvec xi, arma::uvec xp, arma::mat lam);

arma::mat mysum(int n, arma::uvec xi, arma::uvec xp, arma::mat lam);

arma::mat myprod(int n, arma::uvec xi, arma::uvec xp, arma::mat lam);

arma::mat myprod_skip(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, int start, int end);

arma::mat sweep(int N, arma::uvec xi, arma::uvec xp, arma::vec W);

arma::vec sweep2(arma::uvec xi, arma::uvec xp, arma::vec W, arma::vec y);

arma::mat mat_digamma(arma::mat a);

arma::vec vec_digamma(arma::vec a);
