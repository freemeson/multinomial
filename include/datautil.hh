#ifndef DATAUTIL_HH
#define DATAUTIL_HH

struct chi2sigma{
  double chi2;
  double chi2_diff;
  double sigma;
  double sigma_diff;
  unsigned int N_param;
  unsigned int order;
};

struct treestat{
  treestat()
  {orders.clear(); levels.clear(); N_param = 0;}
  void operator+=(treestat rhs){
    if(orders.size()<rhs.orders.size())
      {
	orders.resize(rhs.orders.size(),0);
      }
    for(int i=0; i!= rhs.orders.size(); i++ )
      {
	orders[i]+= rhs.orders[i];
      }


    if(levels.size()<rhs.levels.size())
      {
	levels.resize(rhs.levels.size(),0);
      }
    for(int i=0; i!= rhs.levels.size(); i++ )
      {
	levels[i]+= rhs.levels[i];
      }

    N_param+=rhs.N_param;
  }
  std::vector<unsigned int> orders;
  std::vector<unsigned int> levels;
  unsigned int N_param;
};

template <class T>
struct datapoint{
  datapoint(std::vector<T> const & in_x, T const & in_target, T const & in_weight){
    x = in_x;
    target = in_target;
    weight = in_weight;
  }
  datapoint(T const & in_x, T const & in_target, T const & in_weight){
    x.push_back(in_x);
    target = in_target;
    weight = in_weight;
  }
  ~datapoint(){}
  std::vector<T> x;
  T target;
  T weight;
};

#endif //DATAUTIL_HH
