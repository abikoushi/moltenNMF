void up_A(arma::mat & alpha,
          arma::vec & R,
          const arma::mat & loglambda,
          const arma::vec & y,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const double & a);

void up_A_sp(arma::mat & alpha,
             arma::vec & R,
             const arma::mat & loglambda,
             const arma::vec & yv,
             const arma::uvec & yi,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const double & a);

void up_As_sp(arma::mat & alpha,
              arma::vec & R,
              const arma::mat & loglambda,
              const arma::vec & yv,
              const arma::uvec & yi,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const double & a, 
              const double & NS);

void up_As_sp2(arma::mat & alpha,
               arma::vec & R,
               const arma::mat & loglambda,
               const arma::vec & yv,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const double & a,
               const double & NS);
