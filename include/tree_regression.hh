#ifndef TREE_REGRESSION_HH
#define TREE_REGRESSION_HH
#include <vector>
#include "monopoly.h"
#include "tensor_serie.hh"
#include "full_correlation_tensor_serie.hh"
#include "datautil.hh"

template <class T>
class tree_regression{
  typedef Eigen::LDLT< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > > decomposition;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > Matrix;
public:
  tree_regression(unsigned int dimensions, unsigned int max_order);
  ~tree_regression();

  void fill(std::vector<T> const & x, T const target, T const weight);
  void set(std::vector<datapoint<T> const *> const & data);
  T eval(std::vector<T> const & x);
  T eval(tensor_serie<T> const & x);
  //  void solve(unsigned int f_order, bool solve = false);
  void train();
  void release();
private:
  unsigned int m_max_test_order;
public:
  full_correlation_tensor_serie< T > m_stat;
private:
  tensor_serie<T> m_x;
  bool m_data_owner;
  std::vector<datapoint<T> const * > m_p_data;
  tensor_serie_function<T> const * m_F;  
  tensor_serie_function<T> const * m_pr_axis;
  tensor_serie<T> * m_x_pr_axis_eval;
  unsigned int m_max_order;

  unsigned int m_dimensions;
  bool m_split;
  std::vector<tree_regression<T> *> m_subregion;


};



#endif //TREE_REGRESSION_HH
