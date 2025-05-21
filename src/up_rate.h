void up_B(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const double & b);

void up_B_sp(const int & N,
             arma::mat & beta,
             arma::mat & V,
             arma::mat & logV,
             const arma::mat & alpha,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const arma::uvec & varind,
             const arma::vec & probX0,
             const double & N0,
             const double & b);

void up_B_sp2(const int & N,
              arma::mat & beta,
              arma::mat & V,
              arma::mat & logV,
              const arma::mat & alpha,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const arma::uvec & xp0,
              const double & b);

void up_Bs(const int & N,
           arma::mat & beta,
           arma::mat & V,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & b);


void up_Bs_sp(const int & N,
              arma::mat & beta,
              arma::mat & V,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const arma::vec & probX0,
              const double & N0, const double & NS,
              const double & b);

void up_Bs_sp2(const int & N,
               arma::mat & beta,
               arma::mat & V,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const arma::uvec & xp0,
               const double & N0, const double & NS,
               const double & b);

void up_Bs_sp3(const int & N,
               arma::mat & beta,
               arma::mat & V,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const double & rho,
               const double & N1S,
               const double & NS,
               const double & b);

// void up_B_ww(const int & N,
//              arma::mat & beta,
//              arma::mat & V,
//              arma::mat & logV,
//              const arma::mat & alpha,
//              const arma::uvec & xi,
//              const arma::uvec & xp,
//              const arma::uvec & varind,
//              const arma::vec & weight,
//              const double & b);