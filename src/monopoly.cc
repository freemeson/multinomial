#include "../include/monopoly.h"
#include <iostream>
#include <math.h>
#include <cmath>

#ifdef USE_LAPACK
#include <Accelerate/Accelerate.h>
#endif

//extern tripoly::bool do_log;
//extern tripoly::bool SHOW_LOG;
//extern tripoly::bool SHOW_INFO;

//global variables
//std::vector<double> factorial;

//template <class T> std::vector<std::vector<unsigned long int> > symmetric_tensor<T>::binomial;


binomial_singleton *binomial_singleton::s_instance = 0;

binomial_singleton::binomial_singleton(int max)
{
  //std::cout << "binomial_singleton::binomial_singleton("<< max << ");"<< std::endl;
  init(max);
}

binomial_singleton::~binomial_singleton()
{
  //std::cout << "binomial_singleton::~binomial_singleton();"<< std::endl;
}

double binomial_singleton::double_binomial(int n,int k)
{
	//I've used ROOT's TMath::Binomial here
   if (n<0 || k<0 || n<k) return NAN;
   if (k==0 || n==k) return 1;

   int k1=fmin(k,n-k);
   int k2=n-k1;
   double fact=k2+1;
   for (double i=k1;i>1.;--i)
      fact *= (k2+i)/i;
   return fact;
}

void binomial_singleton::init(int max)
{

  std::vector<unsigned long long int> factorial(max+1, 1);

  unsigned long long int value = 1;
  factorial[0]=1;

  for(int i=1; i!=max+1; i++)
    {
      value*=i;
      factorial[i] = value;
    }
  m_binomial.clear();
  m_binomial.resize(max+1, std::vector<unsigned long int>(max+1, 0));
  bool warning = false;
  for(int i=0; i!=max+1; i++)
    {
      for(int j=0; j!=i+1; j++) // only fill the upper half of the matrix
	{
	  unsigned long int bigger;
	  int smaller;
	  if(j>i-j)
	    {
	      bigger = j; smaller = i-j;
	    }
	  else
	    {
	      bigger = i-j; smaller = j;
	    }
		
	  unsigned long long int k=i;
	  unsigned long long int l=1;
	  value = 1;
	  while(k!=bigger)
	    {
	      value*=k--;	
	      value=value/(l++);
	    }

	  m_binomial[i][j] = value;

	  if( fabs(m_binomial[i][j] - double_binomial(i,j)) > 1.0)
	    {
	      //std::cout << "m_binomial["<< i << "][" << j << "] = " <<  m_binomial[i][j] << ",";
	      warning = true;
	    }
	}
    }
  if(warning ) {std::cout << "Warning in binomial calculation" << std::endl;}
  m_order = max;
}

binomial_singleton * binomial_singleton::instance(int max_order)
{
  //std::cout << "binomial_singleton::instance(" << max_order<< ");" << std::endl;
  if(!s_instance)
    {
      s_instance = new binomial_singleton(max_order);
    }
  else
    {
      s_instance->increase(max_order);
    }
  return s_instance;
}

void binomial_singleton::increase(int max_order)
{
  //std::cout << "binomial_singleton::increase("<< max_order << ");"<< std::endl;
  if(m_order < max_order)
    {
      init(max_order);
    }
  
}

unsigned long int binomial_singleton::get(int up, int down)
{
  //std::cout << "get(" << up << "," << down <<  ") = " << std::endl;//<< m_binomial[up][down] << std::endl;
  return m_binomial[up][down];
}




template <class T> //, unsigned int N, unsigned int K>
void symmetric_tensor<T>::init(unsigned int dimensions, unsigned int order)
{
  m_order = order;
  m_dimensions = dimensions;
  //the number of free parameters of a symmmetric tensor
  //calc_binomial(m_order+dimensions + 1);
  binomial = binomial_singleton::instance(m_order+dimensions+1);


  // if(binomial.size()< m_order+dimensions +2 )
  //   {
  //     std::cout << "++binomial.size() = " << binomial.size() << " needed = " << m_order+dimensions +2 << std::endl;
  //     calc_binomial(m_order+dimensions + 1);
  //   }
  // else{
  //   std::cout << "- binomial.size() = " << binomial.size() << " < needed = " << m_order+dimensions +2 << std::endl;
    // for(int i=0; i!=binomial.size(); i++)
    //   {
    // 	for(int j=0; j!=binomial.size(); j++)
    // 	  {
    // 	    std::cout<< "b[" << i << "]["<< j <<"] = " << binomial[i][j] << std::endl; 
    // 	  }
    //   }
  //  }
  if(dimensions == 1)
    {
      m_size = 1;
    }
  else
    {
      m_size = binomial->get(order+dimensions-1,order);//binomial[order+dimensions-1][order];//TMath::Binomial(order+dimensions -1, order); 
    }
  //std::cout << "m_size = " << m_size << std::endl;
  m_components.clear();
  m_components.resize(m_size,0.0);
  m_entries = 0;
  unsigned long int pos = 0;
  m_diagonal_position_limits.clear();
  m_diagonal_position_limits.resize(m_dimensions,0);

  if(m_order>1)
    {
      for(int i=1; i!=m_dimensions; i++)
	{
	  pos += binomial->get(m_order-2+m_dimensions-i, m_order-2);//binomial[m_order-2 + m_dimensions-i][m_order-2];
	  m_diagonal_position_limits[i] = pos;
	  //std::cout << "i = " << i << " pos = " << pos << std::endl;
	}
    } 
}

template <class T> //, unsigned int N, unsigned int K>
void symmetric_tensor<T>::index_multiplicity(std::vector<unsigned int> const & index, 
					     std::vector<unsigned int> & monomial) const
{
  //no check for vector sizes, this is an internal function
  for(int i=0; i!=m_dimensions; i++)
    {
      monomial[i]=0;
    }
  for(int i=0; i!=m_order; i++)
    {
      monomial[index[i]]++; // the resulting vector might be stored, to ease x^n computation
    }
}

/*
template <class T>
void symmetric_tensor<T>::calc_binomial(int max)
{
  //std::cout << binomial.size() << " " << std::flush;
  std::vector<unsigned long long int> factorial(max+1, 1);
  unsigned long long int value = 1;
  factorial[0]=1;
  for(int i=1; i!=max+1; i++)
    {
      value*=i;
      factorial[i] = value;
    }
  binomial.clear();
  binomial.resize(max+1, std::vector<unsigned long int>(max+1, 0));
  //std::cout << binomial.size() << " " << std::flush;
  for(int i=0; i!=max+1; i++)
    {

      //std::cout << binomial[i ].size() << " " << std::flush;
 	  
      for(int j=0; j!=i+1; j++) // only fill the upper half of the matrix
	{
	  unsigned long int bigger;
	  int smaller;
	  if(j>i-j)
	    {
	      bigger = j; smaller = i-j;
	    }
	  else
	    {
	      bigger = i-j; smaller = j;
	    }
		
	  unsigned long long int k=i;
	  unsigned long long int l=1;
	  value = 1;
	  while(k!=bigger)
	    {
	      value*=k--;	
	      value=value/(l++);
	    }

	  binomial[i][j] = value;

	  //std::cout << "binomial["<< i << "][" << j << "] = " <<  value << std::endl;
	  
	  // while(k!=bigger)
	  //   {
	  //     value*=k--;
	  //   }
	  // value=value/factorial[smaller];
	  // binomial[i][j] = value;

	  //binomial[i][j] = factorial[i]/(factorial[i-j]*factorial[j])
	  if( TMath::Abs((binomial[i][j]-TMath::Binomial(i,j))) > 1.0)
	    {
	      std::cout << "binomial["<< i << "][" << j << "] = " <<  binomial[i][j] << ",";
	    }
	  //std::cout << (binomial[i][j]-TMath::Binomial(i,j)) << " " ;//<< binomial[i][j] << ",";
	}
      //std::cout << std::endl;
    }
}*/

template <class T>
void symmetric_tensor<T>::update_index(std::vector<unsigned int> & index) const
{
  //should check somewhere if index[0] == dimensions-1 already, that's the last index
  unsigned int j=m_order-1;
  while(index[j] == m_dimensions-1) //find least significant index
    {
      j--;
    }
  index[j]++;
  for(int k=j+1; k!=m_order; k++)
    {
      index[k] = index[j];
    }
}

// template <class T>
// void symmetric_tensor<T>::create_monomials()
// {
//   m_monomials.clear();
//   m_monomials.resize(m_size, std::vector<unsigned int>(m_dimensions) );
//   std::vector<unsigned int > index(m_order, 0);
//   std::vector<unsigned int > monomial(m_dimensions);
//   for(int i=0; i!=m_size; i++)
//     {
//       index_multiplicity(index, monomial);
//       m_monomials[i] = monomial;
//       if(i!=m_size - 1) {update_index(index);}
//     }
//   x_powers.clear();
//   x_powers.resize(m_order+1,std::vector<double>(m_dimensions));
// }

// template <class T>
// void symmetric_tensor<T>::fill(std::vector<T> x)
// {
//   for(int i=0; i!= m_order+1; i++)
//     {
//       for(int j=0; j!=m_dimensions; j++)
// 	{
// 	  x_powers[i][j] = TMath::Power(x[j],i);
// 	}
//     }
//   for(int i=0; i!=m_size; i++)
//     {
//       double value = 1.0;
//       for(int j=0; j!=m_dimensions; j++)
// 	{
// 	  if(m_monomials[i][j]!=0)
// 	    {value*= x_powers[m_monomials[i][j]] [j];}
// 	}
//       m_components[i]+=value;
//     }
//   m_entries++;
// }

template <class T>
void symmetric_tensor<T>::create_diad(std::vector<T> const & x, std::vector<T> const & x_diad)
{
  int k=0; 
  for(int i=0; i!=m_dimensions; i++)
    {
      for(int j=m_diagonal_position_limits[i]; j!=x_diad.size(); j++ )
	{
#ifndef LOGDIAD	  
	  m_components[k++] = x[i]*x_diad[ j ] ;
#else
	  m_components[k++] = x[i]+x_diad[ j ] ;
#endif
	  //	  std::cout << "i = " << i << " j = " << j << " k = " << k  << " m_size == "<< m_size << std::endl;
	}
    }
}

template <class T>
void symmetric_tensor<T>::create_diad_1dim(const T & x, const T & x_diad)
{
#ifndef LOGDIAD	  
	  m_components[0] = x*x_diad;
#else
	  m_components[0] = x+x_diad;
#endif
	  //	  std::cout << x << " " << x_diad << std::endl;
}


template <class T>
void symmetric_tensor<T>::fill_diad(std::vector<T> const & x_diad, T const & target, T const & weight)
{
  T value = target*weight;
  if(value == 1)
    {
      for(int i=0; i!=m_size; i++)
	{
#ifndef LOGDIAD
	  m_components[i] += x_diad[i];
#else
	  m_components[i] += exp(x_diad[i]);
#endif
	}
    }
  else{
    if(value == -1)
      {
	for(int i=0; i!=m_size; i++)
	  {
#ifndef LOGDIAD
	  m_components[i] -= x_diad[i];
#else
	  m_components[i] -= exp(x_diad[i]);
#endif
	  }
      }
    else
      {
	for(int i=0; i!=m_size; i++)
	  {

#ifndef LOGDIAD
	  m_components[i] += value*x_diad[i];;
#else
	  m_components[i] += value*exp(x_diad[i]);
#endif
	  }
      }
  }
  m_entries+=weight;  
}

template <class T>
void symmetric_tensor<T>::fill_diad_1dim(const T & x_diad, T const & target, T const & weight)
{
  T value = target*weight;
  if(value == 1)
    {
#ifndef LOGDIAD
	  m_components[0] += x_diad;
#else
	  m_components[0] += exp(x_diad);
#endif
    }
  else{
    if(value == -1)
      {
#ifndef LOGDIAD
	m_components[0] -= x_diad;
#else
	m_components[0] -= exp(x_diad);
#endif
      }
    else
      {
#ifndef LOGDIAD
	  m_components[0] += value*x_diad;
#else
	  m_components[0] += value*exp(x_diad);
#endif
      }
  }
  m_entries+=weight;  
}



template <class T>
T symmetric_tensor<T>::eval_diad(std::vector<T> const & x_diad) const
{
  T value = 0;
  for(int i=0; i!= m_size ; i++)
    {
#ifndef LOGDIAD
      value += m_components[i] * x_diad[i];
#else
      value += m_components[i] * exp(x_diad[i]);
#endif
    }
  return value;
}

template <class T>
T symmetric_tensor<T>::eval_diad_1dim(const T & x_diad) const
{
  T value = 0;
#ifndef LOGDIAD
  value = m_components[0] * x_diad;
#else
  value = m_components[0] * exp(x_diad);
#endif
  return value;
}


template <class T>
void symmetric_tensor<T>::normalize()
{
  for(int i=0; i!=m_size; i++)
    {
      m_components[i]/=(T)m_entries;
    }
}

template <class T>
void symmetric_tensor<T>::normalize(const T factor)
{
  for(int i=0; i!=m_size; i++)
    {
      m_components[i]/=(T)factor;
    }
}

template <class T>
void symmetric_tensor<T>::add(symmetric_tensor<T> const & tensor_b)
{
   for(int i=0; i!=m_size; i++)
     {
       m_components[i] += tensor_b.m_components[i];
     }
   m_entries += tensor_b.m_entries;
}


// template <class T>
// void symmetric_tensor<T>::sum(const symmetric_tensor<T> & left,
// 			      const symmetric_tensor<T> & right)
// {
//   init(left.m_dimensions, left.m_order);
//   for(int i=0; i!=m_size; i++)
//     {
//       m_components[i] = left.m_components[i] + right.m_components[i];
//     }
// }

// template <class T>
// void symmetric_tensor<T>::difference(const symmetric_tensor<T> & left,
// 				     const symmetric_tensor<T> & right)
// {
//   init(left.m_dimensions, left.m_order);
//   for(int i=0; i!=m_size; i++)
//     {
//       m_components[i] = left.m_components[i] - right.m_components[i];
//     }
// }

template <class T> 
unsigned int symmetric_tensor<T>::index_from_monomial(std::vector<unsigned int > const & monomial) const
{
  unsigned int suborder = m_order;
  unsigned int index = m_size;

  for(unsigned int i=0; i!= m_dimensions; i++)
    {      
      for(unsigned int j=0; j!= monomial[i]; j++)
	{
	  index-= binomial->get(suborder+(m_dimensions - i - 2), suborder ); 
	  suborder--;
	}
    }
  index--;
  return index;
}

template <class T>
T const symmetric_tensor<T>::find_monomial(std::vector<unsigned int > const & monomial) const
{
  unsigned int suborder = m_order;
  unsigned int index = m_size;
  //std::cout << "m_dim = " << m_dimensions << " "; 
  //std::cout << "m_order = " << m_order << " " << std::endl;
  // for(unsigned int i=0; i!= m_dimensions; i++)
//     {
//       std::cout << monomial[i] << " ";
//     }
//   std::cout << std::endl;
  for(unsigned int i=0; i!= m_dimensions; i++)
    {      
      for(unsigned int j=0; j!= monomial[i]; j++)
	{
	  /*std::cout << "suborder = " << suborder 
 		    << " i = " << i
 		    << " suborder+(dim-i)-2 = " 
 		    << (suborder+(m_dimensions - i -2))  << std::endl;
	  std::cout << binomial.size() << " " << std::flush;
	  std::cout << binomial[suborder+(m_dimensions - i - 2) ].size() << " " << std::flush;
 	  std::cout << "b[s+d-i-2][s] = " << binomial[ suborder+(m_dimensions - i - 2) ][ suborder ] << std::endl;	  */
	  index-= binomial->get(suborder+(m_dimensions - i - 2), suborder ); 
	    //binomial[ suborder+(m_dimensions - i - 2) ][ suborder ]; 
	  
	  suborder--;
	}
    }
  index--;

  /*  //check if index is really correct
  std::vector<unsigned int> ind(m_order, 0.0);
  std::vector<unsigned int> ind_mon(m_dimensions, 0.0);
  unsigned int k = 0;
  for(unsigned int i=0; i!= m_size; i++)
    {
      index_multiplicity(ind, ind_mon);
      if(ind_mon == monomial) {break;}
      k++;
      if(i!=m_size-1){update_index(ind);     }
    }
  if(k!=index)
    {
      std::cout << k << "!="<< index << " Something went wrong" << std::endl;
      } 
  //  std::cout << "index = " << index << std::endl; */
  return m_components[index];
}

// template <class T>
// T symmetric_tensor<T>::eval(std::vector<T> & x)
// {
//   for(int i=0; i!= m_order+1; i++)
//     {
//       for(int j=0; j!=m_dimensions; j++)
// 	{
// 	  x_powers[i][j] = TMath::Power(x[j],i);
// 	}
//     }
//   T value = 0.0;
//   for(int i=0; i!=m_size; i++)
//     {
//       T xp = 1.0;
//       for(int j=0; j!=m_dimensions; j++)
// 	{	 
// 	  if(m_monomials[i][j]!=0)
// 	    {xp*= x_powers[m_monomials[i][j]] [j];}
// 	}
//       //value+= m_components[i]*(double)m_multiplicities[i]*xp;
//       value+= m_components[i]*xp;
//       //m_components[i]+=value;
//     }
//   return value;
// }

template <class T>
T symmetric_tensor<T>::get_entries()
{
  return m_entries;
}

template <class T>
void symmetric_tensor<T>::print() const
{
  std::cout << "Order  = " << m_order << " Dimension = " << m_dimensions << " Size = " << m_size << std::endl;
  std::vector<unsigned int > index(m_order, 0);
  for(int i=0; i!= m_size; i++)
    {
      for(int j=0; j!= m_order; j++)
	{
	  std::cout << index[j] << " ";
	}
      std::cout << "value = " << m_components[i] << std::endl;
      if(i!= m_size - 1)
	{
	  update_index(index);
	}
    }
}

/*
template <class T>
void symmetric_tensor<T>::release_binomial()
{
  //std::vector<std::vector<unsigned long int > >().swap(binomial);
  }*/

#ifdef USE_LAPACK
template <class T> 
fortran_matrix<T>::fortran_matrix(unsigned int size)
{
  m_size = size;
  m_matrix.clear();
  m_matrix.resize(m_size*m_size);
}

template <class T>
void fortran_matrix<T>::resize(unsigned int size)
{
  m_size = size;
  m_matrix.resize(m_size*m_size);
}

template <class T>
T & fortran_matrix<T>::operator()(unsigned int i, unsigned int j)
{
  return m_matrix[i*m_size + j];
}

template <class T>
T fortran_matrix<T>::operator()(unsigned int i, unsigned int j) const
{
  return m_matrix[i*m_size + j];
}

template <class T>
T * fortran_matrix<T>::get_array()
{
  if(m_size>0)
    {
      return & m_matrix[0];
    }
  else
    {
      std::cout << "null array" << std::endl;
      return NULL;
    }
}

template <class T>
unsigned int fortran_matrix<T>::get_size()
{
  return m_size;
}
#endif

template <class T>
tensor_series<T>::tensor_series(unsigned int dimensions,unsigned int max_order)
{
  m_max_h_order = max_order;
  m_max_g_order = 2*max_order;

  m_h.clear();
  m_g.clear();  
  m_x.clear();
  m_h.resize(m_max_h_order+1);
  m_g.resize(m_max_g_order+1);
  m_x.resize(m_max_g_order+1);
  m_dimensions = dimensions;
  for(int i=0; i!=m_max_h_order+1; i++)
    {
      //std::cout << max_order << std::endl;
      //std::cout << "Init h " << i << std::endl;
      m_h[i].init(dimensions, i);
    }
  for(int i=0; i!=m_max_g_order+1; i++)
    {
      //      std::cout << "Init g,x " << i << std::endl;
      m_g[i].init(dimensions, i);
      m_x[i].init(dimensions, i);
    }
  //  x.clear(); x.resize(dimensions);

}

// template <class T>
// void tensor_series<T>::prepare_multiplicities()
// {
//   for(int i=0; i!=m_max_order +1 ; i++)
//     {
//       m_tensors[i].create_multiplicity_tensor();
//     }
// }

// template <class T>
// void tensor_series<T>::prepare_monomials()
// {
//   for(int i=0; i!=m_max_order +1 ; i++)
//     {
//       m_tensors[i].create_monomials();
//       //std::cout << "mono i=" << i << std::endl;
//     }
// }

template <class T>
void tensor_series<T>::fill(std::vector<T> const & /*in_x*/ x, T const & target, T const & weight)
{ 
  //for(int i=0; i!= m_dimensions ; i++) {x[i] = in_x[i];}
#ifndef LOGDIAD
	    m_x[0].m_components[0] = 1;
#else
	    m_x[0].m_components[0] = 0;
#endif
	    //  m_x[0].m_components[0] = 1;
  for(int i=0; i!=m_max_g_order+1; i++) 
    {
      //m_tensors[i].fill(x);
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_x[i].m_components = x;
#else
	    for(int j=0; j!=m_dimensions ; j++)
	      {
		m_x[i].m_components[j] = log(x[j]);
	      }
#endif
	  }
	  else
	    {
	      m_x[i].create_diad(m_x[1].m_components, m_x[i-1].m_components);	
	    }
	}
      m_g[i].fill_diad(m_x[i].m_components,1, weight );
      if(i<m_max_h_order+1)
	{
	  m_h[i].fill_diad(m_x[i].m_components,target, weight );
	}
    }
}

template <class T>
void tensor_series<T>::fill_1dim(const T & x, T const & target, T const & weight)
{
#ifndef LOGDIAD
	    m_x[0].m_components[0] = 1;
#else
	    m_x[0].m_components[0] = 0;
#endif
  for(int i=0; i!=m_max_g_order+1; i++) 
    {
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_x[i].m_components[0] = x;
#else
	    m_x[i].m_components[0] = log(x);
#endif
	  }
	  else
	    {
	      m_x[i].create_diad_1dim(m_x[1].m_components[0], m_x[i-1].m_components[0]);	
	    }
	}
      m_g[i].fill_diad_1dim(m_x[i].m_components[0],1, weight );
      if(i<m_max_h_order+1)
	{
	  m_h[i].fill_diad_1dim(m_x[i].m_components[0],target, weight );
	}
    }
}


template <class T>
void tensor_series<T>::normalize()
{
  for(int i=0; i!=m_max_h_order+1; i++)
    {
      m_h[i].normalize();      
    }
  for(int i=0; i!=m_max_g_order+1; i++)
    {
      m_g[i].normalize();      
    }
}

template <class T>
void tensor_series<T>::pretrain_add(const tensor_series<T> & series_b)
{
  for(int i=0; i!= m_max_h_order+1; i++)
    {
      m_h[i].add(series_b.m_h[i]);
      //m_h[i].m_components = series_b.m_h[i].m_components;
    }

  for(int i=0; i!= m_max_g_order+1; i++)
    {
      m_g[i].add(series_b.m_g[i]);
    }
}

// template <class T>
// void tensor_series<T>::sum(const tensor_series<T> & left, 
// 			   const tensor_series<T> & right)
// {
//   for(int i=0; i!=m_max_order+1; i++)
//     {
//       m_tensors[i].sum(left.m_tensors[i],right.m_tensors[i]);
//     }
// }

// template <class T>
// void tensor_series<T>::difference(const tensor_series<T> & left,
// 				  const tensor_series<T> & right)
// {
//   for(int i=0; i!=m_max_order+1; i++)
//     {
//       m_tensors[i].difference(left.m_tensors[i],right.m_tensors[i]);
//     }
// }

template <class T>
//TVectorT<double> & tensor_series<T>::vector()
void tensor_series<T>::vector(unsigned int order)
{
  unsigned int size = 0;
  //for(int i=0; i!=m_max_h_order+1; i++)
  for(int i=0; i!=order+1; i++)
    {
      size += m_h[i].m_size;
    }
  //m_vector.ResizeTo(size);
  m_vector.resize(size);
  //m_eigen_vector.resize(size);
  int k = 0;
  std::cout << "size = " << size << std::endl;
  //  for(int i=0; i!=m_max_h_order+1; i++)
  for(int i=0; i!=order+1; i++)
    {
      for(int j=0; j!= m_h[i].m_size; j++)
	{
	  //std::cout << "k= " << k << std::flush;
	  m_vector[k] = m_h[i].m_components[j];
	  //m_eigen_vector[k] = m_h[i].m_components[j];
	  k++;
	}
    }
  std::cout << "end" << std::endl;
  //  return m_vector;
}

template <class T>
//TMatrixTSym<T> tensor_series<T>::matrix(unsigned int order) //, tensor_series<unsigned int> const & multiplicities)
void tensor_series<T>::matrix(unsigned int order) //, tensor_series<unsigned int> const & multiplicities)
{
  //  std::cout << "order = " << order << std::endl;
  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset(order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      offset[i] = size;
      size += m_g[i].m_size;
    }
  //std::cout << "size = " << size << std::endl;
  //TMatrixTSym<T> out(size);
  //double out[size][size];
#ifdef USE_LAPACK
  m_matrix.resize(size);
  //  std::cout << "matrix size = " << m_matrix.get_size() << " " << size << std::endl;
#else
  m_eigen_matrix.resize(size, size);
#endif


  for(int i=0; i!= order+1; i++)
    {
      int v_size = m_g[i].m_size;// vertical size of the chunk under investigation
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_g[j].m_size; // horizontal size of the chunk under investigation
	  //create a h_size*v_size matrix from m_tensor[i+j]
	  std::fill_n(index_v.begin(),i,0);
	  for(int k=0; k!= v_size; k++) 
	    {
	      std::fill_n(index_h.begin(),j,0);
	      for(int l=0; l!= h_size; l++) 
		{
		  for(int m = 0; m!=m_dimensions; m++)
		    {monomial[m] = 0;}
		  for(int m=0; m!=i ;m++ )
		    {
		      //std::cout << index_v[m] << " ";
		      monomial[index_v[m]]++;
		    }
		  //std::cout << std::endl;
		  for(int m=0; m!=j ;m++ )
		    {
		      //std::cout << index_h[m] << " ";
		      monomial[index_h[m]]++;
		    }
		  //std::cout << std::endl << "j+i = "<< (j+i) << std::endl;
		  /*std::cout << "offset[i] = " << offset[i] 
			    << " k = " << k 
			    << "offset[j] = " << offset[j] 
			    << " l = " << l
			    << "offset[i] + k = " << (offset[i]+k) 
			    << "offset[j] + l = " << (offset[j]+l) 
			    << std::endl;*/
		  /*
		  // *multiplicities only along horizontal index
		  //double value = (double)(multiplicities.m_tensors[i].m_components[k]) * m_tensors[i+j].find_monomial(monomial);
		  //the multiplicities should be here, but when I remove it,
		  //it is like encoding it into the unknown F
		  //then... since F.eval() would neeed it, this is a shortcut!
		  //there is no need to use the multiplicities explicitly!!*/
		  double value =  m_g[i+j].find_monomial(monomial);
		  
		  //std::cout << " ["<< offset[i]+k << "][" << offset[j]+l << "] value = " << value;
		  //out(offset[i]+k, offset[j]+l) = value;
#ifdef USE_LAPACK		  
		  m_matrix(offset[i]+k, offset[j]+l) = value; 
#else 
		  m_eigen_matrix(offset[i]+k,offset[j]+l) = value;  
		  m_eigen_matrix(offset[j]+l,offset[i]+k) = value;  
#endif
		  //std::cout << "m_eigen_matrix(offset[i]+k,offset[j]+l) = " << m_eigen_matrix(offset[i]+k,offset[j]+l) << std::endl;
		  //out(offset[j]+l, offset[i]+k) = value;		  
		  //only used for updating the index, it's a const function
		  if(l!=h_size-1){
		    m_g[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_g[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }
  //  return out;
}

template <class T>
//TMatrixTSym<T> tensor_series<T>::matrix(unsigned int order) //, tensor_series<unsigned int> const & multiplicities)
void tensor_series<T>::matrix_1dim(unsigned int order) //, tensor_series<unsigned int> const & multiplicities)
{
  //std::cout << "order = " << order << std::endl;
  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset(order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      offset[i] = size;
      size += m_g[i].m_size;
    }
#ifdef USE_LAPACK
  m_matrix.resize(size);
#else
  m_eigen_matrix.resize(size,size);
#endif

  for(int i=0; i!= order+1; i++)
    {
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  double value =  m_g[i+j].m_components[0]; //every g has only one component
#ifdef USE_LAPACK
	  m_matrix( i , j ) = value; 
#else
	  m_eigen_matrix(i,j) = value;  
	  m_eigen_matrix(j,i) = value;  
#endif
	  //	  std::cout << "[" << i << "][" << j << "] = " << value << std::endl;
	}
    }
  //  return out;
}

template <class T>
//void tensor_series<T>::set_values(TVectorT<double> const & vector)
void tensor_series<T>::set_values()
{
  m_F.clear();
  // m_F.resize(m_max_h_order+1); 
  m_F.resize(m_max_F_order+1);
  unsigned int offset = 0;
  //double const * p = vector.GetMatrixArray();
  //for(int i=0; i!=m_max_h_order+1; i++)
  for(int i=0; i!=m_max_F_order+1; i++)
    {
      m_F[i].init(m_dimensions, i);
      unsigned int size = m_F[i].m_size;      
      //      m_tensors[i].m_components.assign(p+offset,p+offset+size);
      int k=0; 
      for(int j=offset; j!=offset+size; j++)
	{
	  //m_F[i].m_components[k++] = vector[j];//factorial[i];
	  m_F[i].m_components[k++] = m_vector[j];
	}
      offset+=size;      
    }
}

// template <class T>
// void tensor_series<T>::use_multiplicities(tensor_series<unsigned int> & multiplicities)
// {
//   for(int i=0; i!=m_max_order+1; i++)
//     {
//       m_tensors[i].use_multiplicities(multiplicities.m_tensors[i]); 
//       //std::cout << "multi i=" << i << std::endl;
//     }
// }

template <class T>
T tensor_series<T>::eval(std::vector<T> const & x /*in_x*/)
{
  //  for(int i=0; i!= m_dimensions ; i++) {x[i] = in_x[i];}
  T value = 0.0;
  //  for(int i=0; i!= m_max_h_order+1; i++)
  for(int i=0; i!= m_max_F_order+1; i++)
    {
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_x[i].m_components = x;
#else
	    for(int j=0; j!=m_dimensions ; j++)
	      {
		m_x[i].m_components[j] = log(x[j]);
	      }
#endif
	  }
	  else{
	    m_x[i].create_diad(m_x[1].m_components, m_x[i-1].m_components);
	  }
	  value+=m_F[i].eval_diad(m_x[i].m_components);
	}
      else
	{
	  value+=m_F[i].m_components[0];
	}
    }
  return value;
}

template <class T>
T tensor_series<T>::eval_1dim(const T & x)
{
  T value = 0.0;
  //  for(int i=0; i!= m_max_h_order+1; i++)
  for(int i=0; i!= m_max_F_order+1; i++)
    {
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_x[i].m_components[0] = x;
#else
	    m_x[i].m_components[0] = log(x);
#endif
	  }
	  else{
	    m_x[i].create_diad_1dim(m_x[1].m_components[0], m_x[i-1].m_components[0]);
	  }
	  value+=m_F[i].eval_diad_1dim(m_x[i].m_components[0]);
	}
      else
	{
	  value+=m_F[i].m_components[0]; //the 0th order has 1 element only
	}
    }
  return value;
}

template <class T>
void tensor_series<T>::solve()
{
  m_max_F_order = m_max_h_order;//order_condition();//before normalize
  std::cout << "order condition = " << m_max_F_order << std::endl;
  normalize();
  //print("solver");
  vector(m_max_F_order);
  if(m_dimensions == 1)
    {
      //matrix_1dim(m_max_h_order);
      matrix_1dim(m_max_F_order);
    }
  else
    {
      //matrix(m_max_h_order);
      std::cout << "hi" << std::endl;
      matrix(m_max_F_order);
      std::cout << "bye" << std::endl;
    }

  //#define USE_LAPACK
  
#ifdef USE_LAPACK
  __CLPK_integer n = m_vector.size();//h_v_temp.GetNrows(); 
  __CLPK_integer nrhs = 1; //column=1 vectors 
  char uplo = 'L' ; // lower or upper triangular solver, no need to fill the other part of the matrix!
  double * A = m_matrix.get_array();//G.GetMatrixArray();
  __CLPK_integer lda = n;
  double * b = &m_vector[0];//h_v_temp.GetMatrixArray();
  __CLPK_integer ldb = n;
  __CLPK_integer info;

  std::vector<double > A_copy(A, A+n*n);
  std::vector<double > b_copy(b, b+n);

  std::vector<double > A_copy2(A, A+n*n);
  std::vector<double > b_copy2(b, b+n);
  
  dposv_(&uplo, &n, &nrhs, A, &lda, b, &ldb, &info); //positve hermetian solver
  
  if(info == 0)
    {
      std::cout << "Solved!" << std::endl;
      set_values();       
    }
  else
    {
      //std::cout << "Something went wrong! info = " << info << std::endl;
      //for symmetric matrix solver, non-positive definite
      __CLPK_integer ipiv[n];
      double work[n];
      __CLPK_integer lwork = n;
      dsysv_(&uplo, &n, &nrhs, &A_copy[0], &lda, ipiv, &b_copy[0], &ldb, work, &lwork,  &info);
      //m_vector = b_copy;
      //set_values();

      //double x[n]; //I used b_copy instead
  __CLPK_integer ldx = n;

  //char equed = 'N';
  // uplo, n, nrhs, A, lda
  //double AF[n*n]; //output of dsysv = A_copy
  __CLPK_integer ldaf = n;
  //__CLPK_integer ipiv[n];
  //double S[n]; 
  //ldb, b
  //double x[n]; //one row matrix n*ldx , defined above
  //__CLPK_integer ldx = 1;
  /*double rcond;
  double berr[n];
  __CLPK_integer n_err_bnds = 1;
  double err_bnds_norm[n];
  double err_bnds_comp[n];
  __CLPK_integer nparams = 0;
  double params[1];
  double work2[4*n];
  __CLPK_integer iwork[n];
  //info

  dsyrfsx_(&uplo, &equed, &n, &nrhs, &A_copy2[0], &lda, 
	  AF, &ldaf, ipiv, S, &b_copy2[0], &ldb, 
	  &b_copy[0], &ldx, &rcond, berr, 
	  &n_err_bnds, err_bnds_norm, err_bnds_comp, 
	  &nparams, params, work2, iwork, &info);*/ //not in Accelerate.h :(

  double ferr[n];
  double berr[n];
  double work2[3*n];
   __CLPK_integer iwork[n];
   //info;
   dsyrfs_(&uplo, &n, &nrhs, &A_copy2[0], &lda, 
	  &A_copy[0], &ldaf, ipiv, &b_copy2[0], &ldb, 
	    &b_copy[0], &ldx,ferr, berr, work2, iwork, &info);
	    std::cout << "info = " << info << std::endl;
   m_vector = b_copy;
   set_values();

    }

  
//   DSYRFSX	(	CHARACTER 	UPLO,
// CHARACTER 	EQUED,
// INTEGER 	N,
// INTEGER 	NRHS,
// DOUBLE PRECISION, dimension( lda, * ) 	A,
// INTEGER 	LDA,
// DOUBLE PRECISION, dimension( ldaf, * ) 	AF,
// INTEGER 	LDAF,
// INTEGER, dimension( * ) 	IPIV,
// DOUBLE PRECISION, dimension( * ) 	S,
// DOUBLE PRECISION, dimension( ldb, * ) 	B,
// INTEGER 	LDB,
// DOUBLE PRECISION, dimension( ldx, * ) 	X,
// INTEGER 	LDX,
// DOUBLE PRECISION 	RCOND,
// DOUBLE PRECISION, dimension( * ) 	BERR,
// INTEGER 	N_ERR_BNDS,
// DOUBLE PRECISION, dimension( nrhs, * ) 	ERR_BNDS_NORM,
// DOUBLE PRECISION, dimension( nrhs, * ) 	ERR_BNDS_COMP,
// INTEGER 	NPARAMS,
// DOUBLE PRECISION, dimension( * ) 	PARAMS,
// DOUBLE PRECISION, dimension( * ) 	WORK,
// INTEGER, dimension( * ) 	IWORK,
// INTEGER 	INFO 
// )	

#else

  std::cout << "Solving with eigen" << std::endl;
  int n = m_vector.size();//h_v_temp.GetNrows(); 
  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(n);
  for(int i=0; i!=n ; i++)
    {v[i] = m_vector[i];}
  //Eigen::VectorXd x = 
  Eigen::Matrix<T, Eigen::Dynamic, 1> x = 
    //m_eigen_matrix.ldlt().solve(v);
    m_eigen_matrix.partialPivLu().solve(v);
    //m_eigen_matrix.colPivHouseholderQr().solve(v);
    //m_eigen_matrix.fullPivHouseholderQr().solve(v);
    //m_eigen_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(v);
  std::cout << "finished" << std::endl;
  //m_vector.assign(&x[0], &x[n-1]); 
  for(int i=0; i!=n ; i++)
    {
      m_vector[i] = x[i];
      // std::cout << m_vector[i] << std::endl;
    }
  set_values();
#endif
}

template <class T>
unsigned int tensor_series<T>::order_condition()
{
  T n_entries = m_g[0].get_entries();
  std::cout << "n_entries = " << n_entries << std::endl;  
  
      int k = 0;
      double cumulated_variance = 0.0;
      for(int i=1; i< m_max_h_order+1; i++)
	{
	  //for(int j=0; j!=m_dimensions ; j++)
	  //{
	  cumulated_variance += (T)(i*i*i)*m_g[2*i].m_components[0]/((T)n_entries - 1.0)/(T)n_entries;
	  std::cout << "var[" << i << "] = " << sqrt(cumulated_variance) << std::endl;
	      //}
	  if(sqrt(cumulated_variance) < 0.4) //assumes x in [-1, 1]
	    { k++; } 
	  else{break;}
	}
      return k;
  
      //else{return m_max_h_order;}
}

template <class T>
void tensor_series<T>::print(std::string name)
{
  std::cout << "___" << std::endl;
  std::cout << "h Name = '"<< name << "' Max order = " << m_max_h_order << std::endl;
  for(int i=0 ; i!=m_max_h_order+1; i++)
    {
      m_h[i].print();
    }
  std::cout << "---" << std::endl;
  std::cout << "G Name = '"<< name << "' Max order = " << m_max_g_order << std::endl;
  for(int i=0 ; i!=m_max_g_order+1; i++)
    {
      m_g[i].print();
    }
  std::cout << "---" << std::endl;
 
}

template <class T>
void tensor_series<T>::release()
{
  std::vector<symmetric_tensor<T> > empty1(0);
  m_h.swap(empty1);
  std::vector<symmetric_tensor<T> > empty2(0);
  m_g.swap(empty2);
  //m_h.clear();
  //m_g.clear();
}

template class symmetric_tensor<float>;
template class symmetric_tensor<double>;
template class symmetric_tensor<long double>;

template class tensor_series<float>;
template class tensor_series<double>;
template class tensor_series<long double>;
