void up_B(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const double & b);

void up_Bs(const int & N,
           arma::mat & beta,
           arma::mat & V,
           arma::mat & logV,
           const arma::mat & alpha,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & N1,
           const double & b,
           const double & rho);

void up_Bs_sp(const int & N,
               arma::mat & beta,
               arma::mat & V,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const double & rho,
               const double & N1S,
               const double & NS,
               const double & b);

