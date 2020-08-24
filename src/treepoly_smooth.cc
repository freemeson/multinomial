#include "../include/treepoly_smooth.hh"
#include <iostream>
//#include "TMath.h"
#include <cmath>
#include <new>

template<class T>
treepoly_smooth<T>::~treepoly_smooth()
{
  //std::cout << "treepoly_smooth::~treepoly_smooth();" << std::endl;
  //if(m_dimensions!= 1) {delete l3_treepoly_smooth;}
  clear_data();
  if(deep_layers.size()!=0)
    {
      //std::cout << "deleting deep layers" << std::endl;
      delete deep_layers[0];
      delete deep_layers[1];
    }
  std::vector<treepoly_smooth<T> * >().swap(deep_layers);
}

template<class T>
treepoly_smooth<T>::treepoly_smooth(int dimensions, int max_level, int order):
  m_order(order),
  m_shaping_level(5),
  l2(dimensions, order),
  m_x(dimensions, order),
  neural_cut(1e-5),
  do_l4(true),
  sigmoid_width(0.3)
{
  m_dimensions = dimensions;
  m_max_level = max_level;
  avg_t = 0;
  avg_tt = 0;
  sum_w = 0;
  deep_layers.clear();
  skip_me = false;
  per_pi_half = ((T)2.0/(T)3.14159265358979312);
  per_pi = ((T)1.0/(T)3.14159265358979312);
  l3_treepoly_smooth = 0;
  l4 = 0;
  core = 0;
  m_reserved_capacity = 0;
  input = new std::vector<std::vector<T> >;
  n_events = 0;
}

template<class T>
void treepoly_smooth<T>::reserve(int capacity, bool need_input_array)
{
  m_reserved_capacity = capacity;
  try{
    // (*input).reserve(capacity);
    // t.reserve(capacity);
    // w.reserve(capacity);
    if(need_input_array) // in deep_train the input array is provided by the pointer p_input
      {
	if(m_dimensions!=1){
	  (*input).resize(capacity); // i hate to use resize... but these are not complex objects;
	}else{
	  (*input).resize(capacity, std::vector<T>(1));
	}
      }
    t.resize(capacity);
    w.resize(capacity); 
  }
  catch(std::bad_alloc &ba){
    std::cout << "bad alloc at void treepoly_smooth<T>::reserve , ba = " << ba.what() << std::endl;
  }
}

template<class T>
void treepoly_smooth<T>::fill(std::vector<T> const & x, T const & target, T const & weight )
{
  // (*input).push_back(x);
  // t.push_back(target);
  //w.push_back(weight);  

  if(n_events >= m_reserved_capacity) {std::cout << "n_events = " << n_events << " m_reserved_capacity = " << m_reserved_capacity << " input_size = " << input_size << " input.size() = " << (*input).size() << std::endl; }
  (*input)[n_events] = x;
  t[n_events] = target;
  w[n_events] = weight;

  n_events++;
}

template<class T>
void treepoly_smooth<T>::fill_diad_serie(tensor_serie<T> const & x_tensor, T const & target, T const & weight )
{
  // not needed: //  (*input).push_back(x);
  //t.push_back(target);
  //w.push_back(weight);
  if(n_events >= m_reserved_capacity) {std::cout << "n_events = " << n_events << " m_reserved_capacity = " << m_reserved_capacity << " input_size = " << input_size << " input.size() = " << (*input).size() << std::endl; }

  t[n_events] = target;
  w[n_events] = weight;
  n_events++;

  core->fill(x_tensor, 1 ,weight);
  l2.fill(x_tensor, target, weight);
}


template<class T>
void treepoly_smooth<T>::fill_1dim(T const & x, T const & target, T const & weight )
{
  // (*input).push_back(std::vector<T> (1,x));
  // t.push_back(target);
  // w.push_back(weight);
  
  (*input)[n_events][0] = x;
  t[n_events] = target;
  w[n_events] = weight;
  n_events++;
}


template<class T>
void treepoly_smooth<T>::fill(tensor_serie<T> const & x, T const & target, T const & weight )
{
  core->fill(x, 1 ,weight);
  l2.fill(x, target, weight);  
  //avg_t+=target*weight;
  //sum_w+=weight;
}



template<class T>
void treepoly_smooth<T>::set_input_array(std::vector<T> * in_p_input, T * in_p_target, T * in_p_weight, int size)
{
  p_input = in_p_input;
  p_w = in_p_weight;
  p_t = in_p_target;
  input_size = size;
}

template<class T>
void treepoly_smooth<T>::set_input_array(std::vector<T> * in_p_input, int size)
{
  p_input = in_p_input;
  input_size = size;
}

template<class T>
void treepoly_smooth<T>::set_wt_array(T * in_p_target, T * in_p_weight, int size)
{
  p_w = in_p_weight;
  p_t = in_p_target;
  input_size = size;
}

template<class T>
void treepoly_smooth<T>::init_core()
{
  core = new tensor_serie<T>(m_dimensions, 2*m_order);
  //std::cout << "core = " << core << std::endl;
}

template<class T>
void treepoly_smooth<T>::clear_data()
{
  //(*input).clear();w.clear();t.clear();  //this does not touch the pointers
  //std::vector<std::vector<T> >().swap((*input));
  if(input!=0) {delete input;}
  std::vector<T>().swap(w);
  std::vector<T>().swap(t);
}

template<class T>
void treepoly_smooth<T>::train(bool in_posterior_correction)
{
  posterior_correction = in_posterior_correction;
  if(!posterior_correction) {m_shaping_level = 0;do_l4 = false;}
  //the order of the input is not preserved, so better to keep all
  // std::vector<std::vector<T> > input_copy;
  // std::vector<T>  t_copy;
  // std::vector<T>  w_copy;
  // if(m_dimensions != 1)
  //   {
  //     input_copy = input;
  //     t_copy = t;
  //     w_copy = w;
  //   }
  tensor_serie<T> x(m_dimensions, 2*m_order);
  init_core();
  for(int i=0; i!= (*input).size(); i++)
    {
      x.create_diad((*input)[i]);
      fill(x, t[i],w[i]);
    }
  set_input_array(&(*input)[0], (*input).size());
  bool keep_data = m_dimensions != 1 && do_l4 && m_max_level>1;
  solve(x, keep_data,posterior_correction);
  //if(m_dimensions != 1 && do_l4) {l4 = new treepoly_smooth<T>(1,m_shaping_level);}  
  if(posterior_correction && m_dimensions != 1 && do_l4 && m_max_level>1) {
    l4 = new treepoly_smooth<T>(1,4, 5);
    l4->reserve(m_reserved_capacity);
    for(int i=0; i!=(*input).size(); i++)
      {
	m_x.create_diad((*input)[i]);
	l4->fill_1dim(eval(m_x), t[i], w[i]);   
      }
    l4->train();
  }
  //  std::cout<< "end of train(), deleting x? m_dim = "<< m_dimensions << " m_order = " << 2*m_order << std::endl;
  clear_data();

}

template<class T>
void treepoly_smooth<T>::solve(tensor_serie<T> & x, bool keep_data , bool in_posterior_correction)
{
  posterior_correction = in_posterior_correction;
  //set_input_array(&(*input)[0],&t[0],&w[0],(*input).size());
  set_wt_array(&t[0],&w[0],w.size());
  
  std::cout << "dimensions = " << m_dimensions << " level = " << m_max_level << " size = " << input_size << " order = " << m_order << std::endl;
  if(input_size < 1000) {skip_me = true; std::cout << "low statistics, ending the training" << std::endl; return;}
  core->decompose(m_order);
  l2.solve(core->m_decomposed_matrix);
  delete core;
  //  std::cout << "deleted core successfully" << std::endl;
  response.clear();
  response.resize(input_size);

  if(posterior_correction && m_dimensions!=1 && m_shaping_level>0) {l3_treepoly_smooth = new treepoly_smooth(1, m_shaping_level, 3); l3_treepoly_smooth->reserve(m_reserved_capacity);}
  for(int i=0; i!= input_size; i++)
    {
      m_x.create_diad(p_input[i]);
      response[i] = l2.eval(m_x);

      if(posterior_correction && m_dimensions!=1 && m_shaping_level>0 ) {l3_treepoly_smooth->fill_1dim(response[i], p_t[i], p_w[i]);}

    }

  if(posterior_correction && m_dimensions!=1 && m_shaping_level>0) {l3_treepoly_smooth->train();}
  for(int i=0; i!= input_size; i++)
    {
      if(posterior_correction && m_dimensions!=1 && m_shaping_level>0){response[i] = l3_treepoly_smooth->eval_1dim(response[i]);}
      //std::cout << "response[i]" << response[i] << std::endl;
      avg_t +=response[i]*p_w[i];
      avg_tt +=response[i]*response[i]*p_w[i];
      sum_w +=p_w[i];      
    }

  avg_t*=1.0/sum_w; //(T)input_size;//(*input).size();
  avg_tt = avg_tt/sum_w - avg_t*avg_t;
  if(avg_tt == 0.0) {skip_me = true; std::cout << "zero target variance^2, ending training" << std::endl; return;}
  if(avg_tt <= 0.0) {skip_me = true; std::cout << "negative target variance^2, ending training" << std::endl; return;}
  avg_tt = sqrt(avg_tt)/sigmoid_width;//*per_pi_half;
  
  
  if(m_max_level!=1)
    {
      deep_train(x, posterior_correction);
    }
  
  if(!keep_data)
    {
      if(input!=0) {std::vector<std::vector<T> >().swap((*input));}
      std::vector<T>().swap(w);
      std::vector<T>().swap(t);
    }
}

template<class T>
void treepoly_smooth<T>::deep_train(tensor_serie<T> & x, bool in_posterior_correction)
{
  posterior_correction = in_posterior_correction;
  //std::cout << "deep train, avg_t = " << avg_t << " avg_tt = " << avg_tt << std::endl;
  //check max_level here
  //std::cout << "deep train m_order = " << m_order << " m_max_level = " << m_max_level << " m_dimensions = " << m_dimensions << std::endl; 
  //deep_layers.resize(2, treepoly_smooth<T>(m_dimensions,m_max_level-1, m_order) );
  deep_layers.push_back(new treepoly_smooth<T>(m_dimensions,m_max_level-1, m_order));
  deep_layers.push_back(new treepoly_smooth<T>(m_dimensions,m_max_level-1, m_order));
  //std::cout << "deep layers resized" << std::endl;

  //std::cout << "input size = " << input_size << std::endl; 

  deep_layers[0]->init_core();
  deep_layers[0]->reserve(m_reserved_capacity, false);
  deep_layers[1]->init_core();
  deep_layers[1]->reserve(m_reserved_capacity, false);
  int low_index = 0, high_index = input_size-1;
  for(int i=0; i!= input_size; i++)
    {
      x.create_diad(p_input[i]);
      T resp =  neural_response(response[i]);
      T neural_resp = neural_weight(response[i]);
      //      if(resp<avg_t)
      //      std::cout << neural_resp << std::endl;
      if(neural_resp>neural_cut)
	{
	  deep_layers[0]->fill_diad_serie(x,p_t[i]-resp,p_w[i]*neural_resp);	  
	}
      //      else
      if(1.0 - neural_resp>=neural_cut)      
	{
	  deep_layers[1]->fill_diad_serie(x,p_t[i]-resp,p_w[i]*(1.0-neural_resp));
	}
    }
  //response.clear();
  std::vector<T>().swap(response);
  deep_layers[0]->set_input_array(p_input, input_size);
  deep_layers[1]->set_input_array(p_input, input_size);
  if(m_dimensions == 1 ) 
    {//w.clear(), t.clear();
      std::vector<T>().swap(w);std::vector<T>().swap(t);
    }
  deep_layers[0]->solve(x, false, posterior_correction);
  if(posterior_correction && m_dimensions != 1 && m_shaping_level>0) {deep_layers[0]->clear_data();}
  deep_layers[1]->solve(x, false, posterior_correction);
  if(posterior_correction && m_dimensions != 1 && m_shaping_level>0) {deep_layers[1]->clear_data();}
}

template<class T>
T treepoly_smooth<T>::eval(std::vector<T> const & x)
{
  m_x.create_diad(x);
  T response = eval(m_x);

  if(posterior_correction && m_dimensions != 1 && do_l4 && m_max_level>1) {
      return l4->eval_1dim(response);
  }
  else
    {return response;}

}

template<class T>
T treepoly_smooth<T>::eval_1dim(T const & x)
{
  m_x.create_diad_1dim(x);
  T response = eval(m_x);
  return response;
}


template<class T>
T treepoly_smooth<T>::eval(tensor_serie<T> const & x)
{
  if(skip_me) {return avg_t;}
  T value = l2.eval(x);
  if(posterior_correction && m_dimensions!=1 && m_shaping_level>0) {value = l3_treepoly_smooth->eval_1dim(value);}
  //T neural_resp = TMath::ATan( (value-avg_t)/avg_tt);
  
  T neural_resp = neural_weight(value);//TMath::ATan( (value-avg_t)/avg_tt/per_pi_half/sigmoid_width*5.0)*per_pi + 0.5;
  value = neural_response(value);//neural_resp*avg_tt + avg_t;
  //  std::cout << neural_resp << std::endl;
  //std::cout << m_max_level << " , " << value << std::endl;
  if(m_max_level>1)
    {
      //if(value<avg_t)
      if(neural_resp > neural_cut)
	{
	  value += deep_layers[0]->eval(x)*neural_resp; //deep_layers[1]->l2.eval(m_x);
	}
      //else
      if(1.0 - neural_resp >=  neural_cut)
	{
	  value += deep_layers[1]->eval(x)*(1.0-neural_resp); 
	  //value = -0.99;	  
	}
    }
  // else{std::cout << "value 1 = " << value << std::endl;}
  return value;
}

template<class T>
T treepoly_smooth<T>::neural_response(T value)
{
  return atan( (value-avg_t)/avg_tt)*avg_tt+avg_t;
} 

template<class T>
T treepoly_smooth<T>::neural_weight(T value)
{
  return atan( (value-avg_t)/avg_tt/per_pi_half/sigmoid_width*1.0)*per_pi + 0.5;
} 



template<class T>
void treepoly_smooth<T>::pretrain_add(treepoly_smooth<T> const untrained_single_sample)
{
  //Needed for treepoly_smooth_fractions
  //one must take care of normalisations prior to calling this
  //l2.pretrain_add(untrained_single_sample.l2);
  //here should be the addition of the samples, as it is needed for l3 training
  (*input).insert((*input).end(), (*untrained_single_sample.input).begin(), (*untrained_single_sample.input).end());
  w.insert(w.end(), untrained_single_sample.w.begin(), untrained_single_sample.w.end());
  t.insert(t.end(), untrained_single_sample.t.begin(), untrained_single_sample.t.end());
}

template<class T>
T treepoly_smooth<T>::normalize(T norm)
{
  //Needed for treepoly_smooth_fractions
  T wsum = 0;
  for(int i=0; i!= w.size(); i++)
    {
      wsum += w[i];
    }
  for(int i=0; i!= w.size(); i++)
    {
      w[i]*=norm/wsum;
    }
  return wsum;
}

template class treepoly_smooth<double>;
template class treepoly_smooth<long double>;
