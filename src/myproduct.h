arma::mat mysum_t(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::vec myprodvec(const int & n, const arma::uvec & xi, const arma::uvec & xp, const arma::vec & lam);

arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam);

arma::mat myprod(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::vec summyprod(const int & n, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

/*
arma::mat myprod_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam);

arma::mat myprod2(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, arma::vec xx);

arma::mat myprod_skip(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam, const int & start, const int & end);

arma::mat myprod_skip_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam, const int & start, const int & end);

arma::mat myprod_one(const int & N,
                     const arma::uvec & xi,
                     const arma::uvec & xp,
                     const arma::mat & lam,
                     const int & start,
                     const int & end);

arma::vec summyprod2(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, arma::vec xx);

arma::mat mysum_t_r(const int & N, const arma::uvec & xj, const arma::uvec & xp, const arma::mat & lam);

arma::mat mysum(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);
*/

/*
arma::mat sweep(int N, arma::uvec xi, arma::uvec xp, arma::vec W);

arma::vec sweep2(arma::uvec xi, arma::uvec xp, arma::vec W, arma::vec y);
*/
