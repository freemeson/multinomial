#ifndef FULL_CORRELATION_TENSOR_SERIE_HH
#define FULL_CORRELATION_TENSOR_SERIE_HH
#include <vector>
#include <string>
#include "monopoly.h"

template <class T> class tensor_serie;


template <class T>
class tensor_serie_function{
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
public:
  tensor_serie_function(unsigned int dimensions, unsigned int max_order);
  tensor_serie_function(tensor_serie_function<T> const & input);
  ~tensor_serie_function() {};
  T const eval(tensor_serie<T> const & x) const;

  std::vector<symmetric_tensor <T> > m_F;
  

public:
  T m_chi2;
  T m_bias;
  T m_biaso;
  unsigned int m_n_params;
  Matrix m_F_derivate;
  Matrix m_Cov;
  unsigned int m_max_order;
  unsigned int m_dimensions;
};

template <class T>
class full_correlation_tensor_serie{
  friend class full_correlation_tensor_serie<T>;
#ifdef USE_LAPACK
#else
  typedef Eigen::LDLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
#endif

public:
  full_correlation_tensor_serie(unsigned int dimensions,unsigned int max_order);
  ~full_correlation_tensor_serie();
  void normalize();
  void decompose(int order);
  tensor_serie_function<T> const solve(unsigned int const order);//decomposition & ext_decomposed_matrix);
  void create_diad(std::vector<T> const &  x);
  void fill(tensor_serie<T> const & x, T const & target, T const & weight);
  T eval(tensor_serie<T> const & x);
  void create_diad_1dim(const T & x);
  void fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight);
  //T eval_1dim(tensor_serie<T> const & x); //I might have forgotten to define this... seems unimportant
  //  void cov_matrix_1dim(unsigned int order, unsigned int startorder=1); 
  void cov_matrix(unsigned int order); 

  void print(std::string name = "");
  void release();
  tensor_serie_function<T> const * principal_axis();

  std::vector<symmetric_tensor<T> > m_x;
  std::vector<symmetric_tensor<T> > m_xt;
  std::vector<symmetric_tensor<T> > m_Cov_xx;
  std::vector<symmetric_tensor<T> > m_Cov_xxt;
  std::vector<symmetric_tensor<T> > m_Cov_xtxt;


private:
  void vector(unsigned int order);
  //  void cov_vector(unsigned int order,unsigned int startorder=1);
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

#endif //FULL_CORRELATION_TENSOR_SERIE_HH
