#include "../include/tree_regression.hh"
#include "../include/tensor_serie.hh"
#include "../include/full_correlation_tensor_serie.hh"
#include <iostream>

#define __MIN_EVENTS_PER_PARAMETER__ 5
#define __SIGMA_ACCURACY__ 3
#define __ADDITIONAL_TEST_ORDERS__ 3

template <class T>
tree_regression<T>::~tree_regression()
{
  if(m_F) {delete m_F; }
  if(m_pr_axis) {delete m_pr_axis; }
  if(m_subregion.size()!= 0 )
    {
      for(int i=0; i!= m_subregion.size(); i++)
	{
	  delete m_subregion[i];
	}
    }

  release();
}


template <class T>
tree_regression<T>::tree_regression(unsigned int dimensions, unsigned int max_order):
  m_max_test_order(max_order+__ADDITIONAL_TEST_ORDERS__),
  m_stat(dimensions, m_max_test_order),
  m_x(dimensions, 4*m_max_test_order),
  m_subregion(0)
{
  m_max_order = max_order;

  m_dimensions = dimensions;
  m_split = false;
  m_F = 0;
  m_pr_axis = 0;
  m_data_owner = false;

}



template <class T>
void tree_regression<T>::fill(std::vector<T> const & x, T const target, T const weight)
{
  m_x.create_diad(x);
  m_stat.fill(m_x, target, weight);  
  m_p_data.push_back(new datapoint<T>(x, target, weight));
  m_data_owner = true;
}


template <class T>
void tree_regression<T>::set(std::vector<datapoint<T> const *> const & data)
{

  m_p_data = data;
  std::cout << " setting " << m_p_data.size() << std::endl;

  for(int i=0; i!= m_p_data.size(); i++)
    {
      m_x.create_diad(m_p_data[i]->x);
      m_stat.fill(m_x, m_p_data[i]->target, m_p_data[i]->weight);
     
    }    
}

template <class T>
void tree_regression<T>::train()
{

  m_stat.normalize();
  //m_D.normalize();
  
  std::vector<tensor_serie_function<T> > tests;
  int k=0;

  for(int i=0; i<m_max_test_order+1; i=i+1)
    {
      tests.push_back(m_stat.solve(i));
    }

  //determine maximal needed order
  m_max_order = 0;
  T min_expected_chi2 = tests[0].m_chi2 + tests[0].m_bias;
  std::cout << "Echi2[" << 0 << "] = " << min_expected_chi2 << " +- " << sqrt(2)*tests[0].m_bias << std::endl;
  if(tests.size()>0)
    {
      for(int i=1; i < tests.size(); i++)
	{
	  T expected_chi2 = tests[i].m_chi2 + tests[i].m_bias;
	  T chi2_plus_thresold = expected_chi2 + (__SIGMA_ACCURACY__*tests[i].m_bias*sqrt(2));	  
	  if( chi2_plus_thresold < min_expected_chi2 
	      && (double)m_p_data.size()/(double)tests[i].m_n_params > __MIN_EVENTS_PER_PARAMETER__ )
	    {
	      m_max_order = i; min_expected_chi2 = chi2_plus_thresold;
	    }
	  std::cout << "Echi2[" << i << "] = " << expected_chi2 << " +- " << sqrt(2)*tests[i].m_bias << std::endl;
	}

    }
  std::cout << " order = " << m_max_order << std::endl;
  if(m_max_order > m_max_test_order-__ADDITIONAL_TEST_ORDERS__) 
    {
      //split and retrain
      m_split = true;
      m_pr_axis = m_stat.principal_axis();

      if(m_pr_axis==0) {std::cout << "pr_axis is NULL" << std::endl;}
      //      m_D.print("m_D");
      //      m_pr_axis->print("m_pr_axis");
      tensor_serie<T> x_1stDeg(m_dimensions, 1);
      m_subregion.clear();
      m_subregion.resize(2);
      m_subregion[0] = new tree_regression<T>(m_dimensions, m_max_test_order-__ADDITIONAL_TEST_ORDERS__);
      m_subregion[1] = new tree_regression<T>(m_dimensions, m_max_test_order-__ADDITIONAL_TEST_ORDERS__);
      std::vector<datapoint<T> const *> p_data_r0; 
      std::vector<datapoint<T> const *> p_data_r1;
      //here could be a check for dimensions==1 and use x_1stDeg.create_diad_1dim
      for(int i=0; i!= m_p_data.size(); i++)
	{
	  x_1stDeg.create_diad(m_p_data[i]->x);
	  if(m_pr_axis->eval(x_1stDeg)<0)
	    {
	      p_data_r0.push_back(m_p_data[i]);
	    }
	  else
	    {
	      p_data_r1.push_back(m_p_data[i]);
	    }
	  
	}
      m_subregion[0]->set(p_data_r0);
      m_subregion[1]->set(p_data_r1);    
      
      //     if(p_data_r0.size() < 100 ) {exit(0);}
      m_subregion[0]->train();
      m_subregion[1]->train();
    }
  else
    {
      unsigned int f_order = m_max_order; 
      unsigned int offset = 0;
      m_F = new tensor_serie_function<T>(tests[f_order]);
      //      m_F->print("m_F");
    }

  //  pr_axis->print();

  
  release();
  
  //  return stats;
}

template <class T>
void tree_regression<T>::release()
{
  if(m_data_owner)
    {
      for(int i=0; i!= m_p_data.size();i++)
	{
	  delete m_p_data[i];
	}
    }
  m_p_data.clear();
  m_data_owner = false;
}

template <class T>
T tree_regression<T>::eval(std::vector<T> const & x)
{
  m_x.create_diad(x);
   if(!m_split)
    {
      if(m_F==0) {std::cout << "m_F == 0" << std::endl;}
      return m_F->eval(m_x);
    }
  else
    {
      if(m_pr_axis==0) {std::cout << "m_pr_axis == 0" << std::endl;}
      if(m_subregion.size()!=2) {std::cout << "subregion.size()!=2" << std::endl;}
      if(m_pr_axis->eval(m_x)<0)
	{
	  if(m_subregion[0]==0) {std::cout << "m_subregion[0] == 0" << std::endl;}
	  return m_subregion[0]->eval(m_x);
	}
      else
	{
	  if(m_subregion[1]==0) {std::cout << "m_subregion[1] == 0" << std::endl;}
	  return m_subregion[1]->eval(m_x);
	}
    }
}

template <class T>
T tree_regression<T>::eval(tensor_serie<T> const & x)
{
  //  std::cout << "eval " << x.m_S[1].m_components[0] << std::endl;
   if(!m_split)
    {
      if(m_F==0) {std::cout << "m_F == 0" << std::endl;}
      return m_F->eval(x);
    }
  else
    {
      if(m_pr_axis==0) {std::cout << "m_pr_axis == 0" << std::endl;}
      if(m_subregion.size()!=2) {std::cout << "subregion.size()!=2" << std::endl;}
      if(m_pr_axis->eval(x)<0)
	{
	  if(m_subregion[0]==0) {std::cout << "m_subregion[0] == 0" << std::endl;}
	  return m_subregion[0]->eval(x);
	}
      else
	{
	  if(m_subregion[1]==0) {std::cout << "m_subregion[0] == 1" << std::endl;}
	  return m_subregion[1]->eval(x);
	}
    }
}



template class tree_regression<float>;
template class tree_regression<double>;
template class tree_regression<long double>;


