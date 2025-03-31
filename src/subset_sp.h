void subset_spx(arma::uvec & xi2, arma::uvec & xp2,
                const arma::uvec & xi, const arma::uvec & xp,
                const arma::uvec & uid);

arma::uvec subset_xp(const arma::uvec & xi, const arma::uvec & xp,
                     const arma::uvec & uid);

void filter_y(arma::vec & yv2, arma::uvec & yi2,
              const arma::vec & yv, const arma::uvec & yi,
              const arma::uvec & ind);
