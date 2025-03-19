
double sumouter2vec(arma::vec & x, arma::vec & y);

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l);

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l,
                    const arma::field<arma::vec> & weight);

double sumouterprod_s(const arma::field<arma::mat> & V,
                      const int & k,
                      const int & l,
                      const arma::field<arma::uvec> & uid);

double sumouterprod_s(const arma::field<arma::mat> & V,
                      const int & k,
                      const int & l,
                      const arma::field<arma::vec> & weight,
                      const arma::field<arma::uvec> & uid);

double sumouterprod_mat(const arma::mat & V,
                        const arma::uvec & varind,
                        const int & k,
                        const int & l,
                        const arma::field<arma::vec> & weight);

double sumouterprod_mat(const arma::mat & V,
                        const arma::uvec varind,
                        const int & k,
                        const int & l,
                        const arma::field<arma::vec> & prob,
                        const arma::vec & freq);