arma::mat mysum_t(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::vec myprodvec(const int & n, const arma::uvec & xi, const arma::uvec & xp, const arma::vec & lam);

arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam);

arma::mat myprod(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::vec summyprod(const int & n, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::mat mysum_t_rv(const int & N,
                     const arma::uvec & xp,
                     const arma::mat & lam);

// arma::mat mysum_t_ww(const int & N,
//                      const arma::uvec & xi,
//                      const arma::uvec & xp,
//                      const arma::mat & lam,
//                      const arma::vec & weight);