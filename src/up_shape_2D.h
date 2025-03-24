double up_A_2D(arma::field<arma::mat> & alpha,
            const arma::field<arma::mat> & logV,
            const arma::vec & y,
            const arma::umat & X,
            const double & a,
            const int & L,
            const arma::uvec & dims,
            const int & k);

double up_As_2D(arma::field<arma::mat> & alpha,
                const arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const double & a,
                const int & L,
                const arma::field<arma::uvec> & uid,
                const int & k);

double up_As_2D(arma::field<arma::mat> & alpha,
                const arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const double & a,
                const int & L,
                const arma::field<arma::uvec> & uid,
                const int & k,
                const double NS);