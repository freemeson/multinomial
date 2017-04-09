#ifndef MONOPOLY_H
#define MONOPOLY_H
#include <vector>
#include <string>

#ifndef USE_LAPACK
#include <Eigen/Dense>
#endif


//global variables
//extern std::vector<double> factorial;
//extern std::vector<std::vector<unsigned long int> > binomial;

/* namespace tripoly */
/* { */
/*   bool do_log = false; */
/*   bool SHOW_LOG=true; */
/*   bool SHOW_INFO=true; */
/* #define tripoly::COUT if(tripoly::do_log){ std::cout  */
/* #define tripoly::LOG if(tripoly::SHOW_LOG) { std::cout  */
/* #define tripoly::INFO if(tripoly::SHOW_INFO) { std::cout  */
/* #define tripoly::ENDL std::endl;} */
/* #define tripoly::FLUSH std::flush;} */
/* } */

//#define LOGDIAD
//#undef LOGDIAD

//class declarations

class binomial_singleton{
private:
  binomial_singleton(int max);
  ~binomial_singleton();
  void init(int max);
  static binomial_singleton * s_instance;
  std::vector<std::vector<unsigned long int > > m_binomial;
  int m_order;
public:
  void increase(int max_order);
  static binomial_singleton * instance(int max_order);
  unsigned long int get(int up, int down);
};

template <class T> //, unsigned int N, unsigned int K>
class symmetric_tensor{
  friend class symmetric_tensor<T>;
public:
  symmetric_tensor(){}
  ~symmetric_tensor(){}
  void init(unsigned int dimensions, unsigned int order);
  //void calc_binomial(int max);
  void update_index(std::vector<unsigned int> & index) const;
  //void create_monomials();
  //void fill(std::vector<T> x);
  void normalize();
  void normalize(const T factor);
//   void sum(const symmetric_tensor<T> & left, const symmetric_tensor<T> & right);
  void add(const symmetric_tensor<T> & tensor_b);
//   void difference(const symmetric_tensor<T> & left, const symmetric_tensor<T> & right);
  unsigned int index_from_monomial(std::vector<unsigned int > const & monomial) const;
  T const find_monomial(std::vector<unsigned int > const & monomial) const;
  //  T eval(std::vector<T> & x);
  void create_diad(std::vector<T> const &x, std::vector<T> const & x_diad);
  void fill_diad(std::vector<T> const & x_diad, T const & target, T const & weight);
  T eval_diad(std::vector<T> const & x_diad) const;

  void create_diad_1dim(const T &x, const T & x_diad);
  void fill_diad_1dim(const T & x_diad, T const & target, T const & weight);
  T eval_diad_1dim(const T & x_diad) const;
  T get_entries();
  //void release_binomial();

  void print() const;
 
  std::vector<T> m_components;
  //std::vector<std::vector<unsigned int> > m_monomials;
  unsigned int m_order;
  unsigned int m_dimensions;
  unsigned int m_size;
  std::vector<unsigned long int> m_diagonal_position_limits; //diag position+1 of lower degree tensors!, so it can be used to limit cycles 

private:
  void index_multiplicity(std::vector<unsigned int> const & index, 
			  std::vector<unsigned int> & monomial) const;
  T m_entries; //because it might be weighted 
  //std::vector<std::vector<double> > x_powers;  
  //std::vector<std::vector<unsigned long int> > binomial;
  binomial_singleton * binomial; 
};


#ifdef USE_LAPACK
template <class T>
class fortran_matrix{
public:
  fortran_matrix(unsigned int size = 0);
  ~fortran_matrix(){};
  void resize(unsigned int size);
  T & operator()(unsigned int i, unsigned int j);
  T operator()(unsigned int i, unsigned int j) const;
  T * get_array();
  unsigned int get_size();
private:
  unsigned int m_size;
  std::vector<T> m_matrix;
};
#endif

template <class T>
class tensor_series{
public:
  tensor_series(unsigned int dimensions,unsigned int max_order);
  ~tensor_series(){}
  //void prepare_multiplicities();
  //void prepare_monomials();
  void fill(std::vector<T> const & x, T const & target, T const & weight);
  void normalize();
  //void sum(const tensor_series<T> & left, const tensor_series<T> & right);
  void pretrain_add(tensor_series<T> const & series_b);
  //void difference(const tensor_series<T> & left, const tensor_series<T> & right);
  //  TVectorT<double> & vector();
  //TMatrixTSym<T> matrix(unsigned int order); //,tensor_series<unsigned int> const & multiplicities); //there is no need for multiplicity tensor anymore  
  //void set_values(TVectorT<double> const & vector);
  void solve();
  //void use_multiplicities(tensor_series<unsigned int> & multiplicities);
  T eval(std::vector<T> const & x);

  void fill_1dim(const T & x, T const & target, T const & weight);
  T eval_1dim(const T & x);

  void print(std::string name = "");
  void release();

  std::vector<symmetric_tensor<T> > m_h;
  std::vector<symmetric_tensor<T> > m_F;
  std::vector<symmetric_tensor<T> > m_g;
  std::vector<symmetric_tensor<T> > m_x;
private:
  void vector(unsigned int order);
  void matrix(unsigned int order); 
  void matrix_1dim(unsigned int order); 
  void set_values();
  unsigned int order_condition();
  unsigned int m_max_h_order;
  unsigned int m_max_g_order;
  unsigned int m_max_F_order;
  unsigned int m_dimensions;
  //TVectorT<double> m_vector;
#ifdef USE_LAPACK
  std::vector<double > m_vector;
  fortran_matrix<double > m_matrix;
#else
  std::vector<T> m_vector;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > m_eigen_matrix;
#endif
  //  std::vector<T> x;
};


#endif //MONOPOLY_H
