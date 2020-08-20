#ifndef TREEPOLY_SMOOTH_HH
#define TREEPOLY_SMOOTH_HH
#include "tensor_serie.hh"

template<class T>
class treepoly_smooth{
  typedef tensor_serie<T> learning_node;
  typedef tensor_series<T> shaping_node;
public:
  treepoly_smooth(int dimensions, int max_level, int order = 2);
  ~treepoly_smooth();
  void init_core();
  void set_input_array(std::vector<T> * p_input, T * p_target, T * p_weight, int size);
  void set_input_array(std::vector<T> * p_input, int size);
  void set_wt_array(T * p_target, T * p_weight, int size);
  void fill(tensor_serie<T> const & x, T const & target, T const & weight );
  void fill(std::vector<T> const & x, T const & target, T const & weight );
  void fill_1dim(T const & x, T const & target, T const & weight );
  void fill_diad_serie(tensor_serie<T> const & x_tensor, T const & target, T const & weight);
  void reserve(int capacity, bool need_input_array = true);
  T neural_response(T resp);
  T neural_weight(T resp);
  T eval(tensor_serie<T> const & x);
  T eval(std::vector<T> const & x);
  T eval_1dim(T const & x);
  void train(bool in_posterior_correction = true);
  void solve(tensor_serie<T> & x, bool keep_data = false, bool in_posterior_correction = true);
  void deep_train(tensor_serie<T> & x, bool in_posterior_correction = true);
  void clear_data();  
  
  int m_dimensions;
  int m_max_level;
  int m_order;
  int m_shaping_level;
  int m_reserved_capacity;
  
  learning_node l2;
  //  shaping_node  l3;
  //shaping_node  * l4;
  treepoly_smooth<T>  * l4;
  treepoly_smooth<T> * l3_treepoly_smooth;
  std::vector<treepoly_smooth<T> * > deep_layers;
  tensor_serie<T> * core;
  tensor_serie<T> m_x;
  T avg_t;
  T avg_tt;
  T sum_w;

  std::vector<std::vector<T> > * input;
  std::vector<T>  t;
  std::vector<T>  w;
  std::vector<T>  response;
  unsigned int n_events; //after the filling this becomes input_size

  std::vector<T> * p_input;
  T  *p_t;
  T  *p_w;
  int input_size;
  bool skip_me;
  T  per_pi_half;
  T  per_pi;
  T neural_cut;
  bool do_l4;
  bool posterior_correction;
  T sigmoid_width;

  void pretrain_add(treepoly_smooth<T> const untrained_single_sample);
  T normalize(T norm = 1.0);
};


#endif //TREEPOLY_SMOOTH_HH
