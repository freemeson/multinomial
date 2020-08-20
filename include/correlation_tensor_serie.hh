#ifndef CORRELATION_TENSOR_SERIE_HH
#define CORRELATION_TENSOR_SERIE_HH
#include <vector>
#include <string>
#include "monopoly.h"

template <class T> class tensor_serie;

template <class T>
class correlation_tensor_serie{
  friend class correlation_tensor_serie<T>;
#ifdef USE_LAPACK
#else
  typedef Eigen::LDLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
#endif

public:
  correlation_tensor_serie(unsigned int dimensions,unsigned int max_order);
  ~correlation_tensor_serie();
  void normalize();
  void pretrain_add(correlation_tensor_serie<T> const & series_b);
  void decompose(int order);
  void solve(decomposition & ext_decomposed_matrix);
  void create_diad(std::vector<T> const &  x);
  void fill(tensor_serie<T> const & x, T const & target, T const & weight);
  T eval(tensor_serie<T> const & x);
  void create_diad_1dim(const T & x);
  void fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight);
  T eval_1dim(tensor_serie<T> const & x);
  void cov_matrix_1dim(unsigned int order, unsigned int startorder=1); 
  void cov_matrix(unsigned int order, unsigned int startorder=1); 
  T sigma_distance(correlation_tensor_serie<T> const & ext_serie);
  void variable_selection(correlation_tensor_serie<T> const & ext_serie, std::vector<T> & sigmas, std::vector<int> & order);

  /*
   //this function is to fit the amplitude/cross section of samples from theoretical distributions
    // to the sample of the data
  void fractions(correlation_tensor_serie<T> const & A, 
		 correlation_tensor_serie<T> const & B, 
		 T & fraction, T & error);
	*/
  void print(std::string name = "");
  void release();

  std::vector<symmetric_tensor<T> > m_S;
  std::vector<symmetric_tensor<T> > m_C; // serie of the correlations

private:
  void vector(unsigned int order);
  void cov_vector(unsigned int order,unsigned int startorder=1);
  void matrix(unsigned int order); 
  void matrix_1dim(unsigned int order); 
  void set_values();
  unsigned int m_max_order;
  unsigned int m_dimensions;
  //TVectorT<double> m_vector;
  std::vector<T> m_vector;
  T m_sum_weight;
  unsigned long int m_entries;
  unsigned int max_significant_order;
#ifdef USE_LAPACK
  fortran_matrix<double > m_matrix;
#else

public:
  Matrix m_eigen_matrix;
  decomposition m_decomposed_matrix;
#endif
  //  std::vector<T> x;
};

#endif //CORRELATION_TENSOR_SERIE_HH
