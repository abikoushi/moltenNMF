arma::mat mysum_t(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);// [[Rcpp::export]]

arma::mat mysum_t_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam);

arma::mat mysum(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::mat myprod(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::mat myprod_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam);

arma::mat myprod2(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, arma::vec xx);

arma::mat myprod_skip(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam, const int & start, const int & end);

arma::mat myprod_skip_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam, const int & start, const int & end);

arma::vec summyprod2(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, arma::vec xx);

arma::mat sweep(int N, arma::uvec xi, arma::uvec xp, arma::vec W);

arma::vec sweep2(arma::uvec xi, arma::uvec xp, arma::vec W, arma::vec y);

arma::mat mat_digamma(arma::mat a);

arma::vec vec_digamma(arma::vec a);
