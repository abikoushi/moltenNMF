
double sumouter2vec(arma::vec & x, arma::vec & y);

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l);

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l,
                    const arma::field<arma::vec> & weight);
