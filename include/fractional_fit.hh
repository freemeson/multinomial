#ifndef FRACTIONAL_FIT_HH
#define FRACTIONAL_FIT_HH
#include <Eigen/Dense>
#include <vector>

/* this class intends to fit probability of processes to the data
 * with sim_1(x), sim_2(x), sim_3(x) samples, a,b,c probabilities, where a+b+c=1
 * and data(x) sample
 *   a*sim_1 + b*sim_2 + c*sim_3 = data
 * is achieved, via fitting the moments of the 1,2,3 distributions to the moments of the data
 * utilized by correlation_tensor_serie
 */


template<class T>
class fractional_fit{
typedef Eigen::LDLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
 public:
 
  fractional_fit(Matrix const & A, Matrix const & B, Matrix const & D, 
		 //Vector const & m, Vector const & n, Vector const & d );
		 std::vector<T> const &m, std::vector<T> const &n, std::vector<T> const &d);
  ~fractional_fit(){}
//  void solve();
  void analytic_solve();
  double operator()(double a);

  T a_min, chi2; 
private:
  Matrix A, B, D;
  Vector m, n, d;
  
  Matrix Sigma; //depends on the 'a' fraction
  Vector l;


};



#endif //FRACTIONAL_FIT_HH
