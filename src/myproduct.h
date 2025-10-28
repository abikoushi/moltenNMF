arma::vec myprodvec(const arma::uword & n, const arma::uvec & xi, const arma::uvec & xp, const arma::vec & lam);

arma::vec myprodvec_sub(const arma::uword & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam);

arma::mat myprod(const arma::uword & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::vec summyprod(const arma::uword & n, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::mat mysum_t(const arma::uword & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam);

arma::mat mysum_t_vec(const arma::uword & N,
                      const arma::uvec & xi,
                      const arma::uvec & xp,
                      const arma::vec & lam);

// arma::mat mysum_t_rv(const int & N,
//                      const arma::uvec & xp0,
//                      const arma::vec & fl);

// arma::mat mysum_t_ww(const int & N,
//                      const arma::uvec & xi,
//                      const arma::uvec & xp,
//                      const arma::mat & lam,
//                      const arma::vec & weight);