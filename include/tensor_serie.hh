#ifndef TENSOR_SERIE_HH
#define TENSOR_SERIE_HH
#include <vector>
#include <string>
#include "monopoly.h"


template <class T>
class tensor_serie{
  friend class tensor_serie<T>;
#ifdef USE_LAPACK
#else
  typedef Eigen::LDLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  //typedef Eigen::FullPivLU< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
#endif

public:
  tensor_serie(unsigned int dimensions,unsigned int max_order);
  ~tensor_serie();
  void normalize();
  void pretrain_add(tensor_serie<T> const & series_b);
  void decompose(int order);
  void solve(decomposition & ext_decomposed_matrix);
  void create_diad(std::vector<T> const &  x);
  void fill(tensor_serie<T> const & x, T const & target, T const & weight);
  T eval(tensor_serie<T> const & x) const;
  void create_diad_1dim(const T & x);
  void fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight);
  T eval_1dim(tensor_serie<T> const & x) const;

  void print(std::string name = "") const;
  void release();
  tensor_serie<T> const * principal_axis();

  std::vector<symmetric_tensor<T> > m_S;
private:
  void vector(unsigned int order);
  void matrix(unsigned int order); 
  void matrix_1dim(unsigned int order); 
  void set_values();
  unsigned int m_max_order;
  unsigned int m_dimensions;
  //TVectorT<double> m_vector;
  std::vector<T> m_vector;
#ifdef USE_LAPACK
  fortran_matrix<double > m_matrix;
#else
 
public:
  Matrix m_eigen_matrix;
  decomposition m_decomposed_matrix;
#endif

  //public:
  //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  //  std::vector<T> x;
};

#endif //TENSOR_SERIE_HH
