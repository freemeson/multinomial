#include "../include/full_correlation_tensor_serie.hh"
#include "../include/tensor_serie.hh"
#include <iostream>
//#include "../include/kissing_ellipsoid_solver.hh"
//#include "TMath.h"


template <class T>
full_correlation_tensor_serie<T>::~full_correlation_tensor_serie()
{
  //  std::cout << "full_correlation_tensor_serie::~full_correlation_tensor_serie() m_max_order = "<< m_max_order << std::endl;
}

template <class T>
full_correlation_tensor_serie<T>::full_correlation_tensor_serie(unsigned int dimensions,unsigned int max_order)
{
  m_max_order = max_order;

  m_x.clear();
  m_x.resize(2*m_max_order+1);
  m_xt.clear();
  m_xt.resize(m_max_order+1);
  m_Cov_xx.clear();
  m_Cov_xx.resize(4*m_max_order+1);
  m_Cov_xxt.clear();
  m_Cov_xxt.resize(3*m_max_order+1);
  m_Cov_xtxt.clear();
  m_Cov_xtxt.resize(2*m_max_order+1);



  m_dimensions = dimensions;
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_xt[i].init(m_dimensions, i);
    }


  for(int i=0; i!=2*m_max_order+1; i++)
    {
      m_x[i].init(m_dimensions, i);
      m_Cov_xtxt[i].init(m_dimensions, i);
    }
  for(int i=0; i!=3*m_max_order+1; i++)
    {
      m_Cov_xxt[i].init(m_dimensions, i);
    }

  for(int i=0; i!=4*m_max_order+1; i++)
    {
      m_Cov_xx[i].init(m_dimensions, i);
    }
  
  m_sum_weight = 0.0;
  m_entries=0;
}

template <class T>
void full_correlation_tensor_serie<T>::fill(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_xt[i].fill_diad(x.m_S[i].m_components,target, weight );      
    }

  for(int i=0; i!=2*m_max_order+1; i++) 
    {
      m_x[i].fill_diad(x.m_S[i].m_components,1., weight );      
      m_Cov_xtxt[i].fill_diad(x.m_S[i].m_components,target*target, weight*weight );
    }

  for(int i=0; i!=3*m_max_order+1; i++) 
    {
      m_Cov_xxt[i].fill_diad(x.m_S[i].m_components,target, weight*weight );
    }

  for(int i=0; i!=4*m_max_order+1; i++) 
    {
      m_Cov_xx[i].fill_diad(x.m_S[i].m_components,1., weight*weight );
    }

  m_sum_weight+=weight;
  m_entries++;
}


template <class T>
void full_correlation_tensor_serie<T>::fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_xt[i].fill_diad_1dim(x.m_S[i].m_components[0],target, weight );      
    }

  for(int i=0; i!=2*m_max_order+1; i++) 
    {
      m_x[i].fill_diad_1dim(x.m_S[i].m_components[0],1., weight );      
      m_Cov_xtxt[i].fill_diad_1dim(x.m_S[i].m_components[0],target*target, weight*weight );
    }

  for(int i=0; i!=3*m_max_order+1; i++) 
    {
      m_Cov_xxt[i].fill_diad_1dim(x.m_S[i].m_components[0],target, weight*weight );
    }

  for(int i=0; i!=4*m_max_order+1; i++) 
    {
      m_Cov_xx[i].fill_diad_1dim(x.m_S[i].m_components[0],1., weight*weight );
    }


  m_sum_weight+=weight;
  m_entries++;
}


template <class T>
void full_correlation_tensor_serie<T>::normalize()
{
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_xt[i].normalize((T)m_entries);
    }
  for(int i=0; i!=2*m_max_order+1; i++)
    {
      m_x[i].normalize((T)m_entries);
      m_Cov_xtxt[i].normalize((T)m_entries);
    }
  for(int i=0; i!=3*m_max_order+1; i++)
    {
      m_Cov_xxt[i].normalize((T)m_entries);
    }
  for(int i=0; i!=4*m_max_order+1; i++)
    {
      m_Cov_xx[i].normalize((T)m_entries);
    }

  //std::cout << "sum of weights = " << m_sum_weight << " entries = " << m_entries << std::endl;
}

template <class T>
void full_correlation_tensor_serie<T>::vector(unsigned int order)
{
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      size += m_xt[i].m_size;
    }
  m_vector.resize(size);
  //m_eigen_vector.resize(size);
  int k = 0;
  //std::cout << "size = " << size << std::endl;  
  for(int i=0; i!=order+1; i++)
    {
      for(int j=0; j!= m_xt[i].m_size; j++)
	{
	  m_vector[k] = m_xt[i].m_components[j];
	  k++;
	}
    }
  //std::cout << "end" << std::endl;
}

/*
template <class T>
void full_correlation_tensor_serie<T>::cov_vector(unsigned int order, unsigned int startorder)
{
  unsigned int size = 0;
  for(int i=startorder; i!=order+1; i++) //this vector starts from i=1 
    {
      size += m_xt[i].m_size;
    }
  m_vector.resize(size);
  //m_eigen_vector.resize(size);
  int k = 0;
  //std::cout << "size = " << size << std::endl;  
  for(int i=startorder; i!=order+1; i++)
    {
      for(int j=0; j!= m_xt[i].m_size; j++)
	{
	  m_vector[k] = m_xt[i].m_components[j];
	  k++;
	}
    }
  //std::cout << "end" << std::endl;
}
*/

template <class T>
void full_correlation_tensor_serie<T>::matrix(unsigned int order)
{  
  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset(order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      offset[i] = size;
      size += m_x[i].m_size;
    }
  
  m_eigen_matrix.resize(size, size);


  for(int i=0; i!= order+1; i++)
    {
      int v_size = m_x[i].m_size;// vertical size of the chunk under investigation
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_x[j].m_size; // horizontal size of the chunk under investigation
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
		  double value =  m_x[i+j].find_monomial(monomial);
		  
		  //std::cout << " ["<< offset[i]+k << "][" << offset[j]+l << "] value = " << value;
		  //out(offset[i]+k, offset[j]+l) = value;

		  m_eigen_matrix(offset[i]+k,offset[j]+l) = value;  
		  m_eigen_matrix(offset[j]+l,offset[i]+k) = value;  

		  //std::cout << "m_eigen_matrix(offset[i]+k,offset[j]+l) = " << m_eigen_matrix(offset[i]+k,offset[j]+l) << std::endl;
		  //out(offset[j]+l, offset[i]+k) = value;		  
		  //only used for updating the index, it's a const function
		  if(l!=h_size-1){
		    m_x[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_x[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }
  //  return out;
}

template <class T>
void full_correlation_tensor_serie<T>::matrix_1dim(unsigned int order)
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
      size += m_x[i].m_size;
    }

  m_eigen_matrix.resize(size,size);

  for(int i=0; i!= order+1; i++)
    {
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  double value =  m_x[i+j].m_components[0]; //every g has only one component

	  m_eigen_matrix(i,j) = value;  
	  m_eigen_matrix(j,i) = value;  
	  //	  std::cout << "[" << i << "][" << j << "] = " << value << std::endl;
	}
    }
  //  return out;
}

template <class T>
void full_correlation_tensor_serie<T>::cov_matrix(unsigned int order)
{  
  //for density tensor series the first order is 1, because the normalization is not a free parameter, therefore it does not show up in the covariance matrix
  //cov_vector(order);

  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset_xt(order+1,0);
  std::vector<unsigned int> offset_xx(2*order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size_xt = 0;
  unsigned int size_xx = 0;
  for(int i=0; i!=order+1; i++) 
    {
      offset_xt[i] = size_xt;
      size_xt += m_Cov_xtxt[i].m_size;
    }

  for(int i=1; i!=2*order+1; i++)
    {
      offset_xx[i] = size_xx;
      size_xx += m_Cov_xx[i].m_size;
    }

  m_eigen_matrix.resize(size_xt + size_xx, size_xt + size_xx);
  m_eigen_matrix.setZero();

  for(int i=0; i!= order+1; i++)
    {
      int v_size = m_Cov_xtxt[i].m_size;// vertical size of the chunk under investigation
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_Cov_xtxt[j].m_size; // horizontal size of the chunk under investigation
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

		  std::cout << m_Cov_xtxt[i+j].find_monomial(monomial) << 
		    " + " << m_xt[i].m_components[k] * m_xt[j].m_components[l] << 
		    " - " << m_Cov_xxt[i].m_components[k] * m_xt[j].m_components[l] << 
		    " - " << m_xt[i].m_components[k] * m_Cov_xxt[j].m_components[l] << std::endl;

		  double value 
		    =  m_Cov_xtxt[i+j].find_monomial(monomial) 
		    +  m_xt[i].m_components[k] * m_xt[j].m_components[l]
		    -  m_Cov_xxt[i].m_components[k] * m_xt[j].m_components[l]
		    -  m_xt[i].m_components[k] * m_Cov_xxt[j].m_components[l];

		  m_eigen_matrix(offset_xt[i]+k,offset_xt[j]+l) = value;  
		  m_eigen_matrix(offset_xt[j]+l,offset_xt[i]+k) = value;  

		  if(l!=h_size-1){
		    m_Cov_xtxt[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_Cov_xtxt[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }


  index_h.resize(2*order+1);
  index_h.clear(); //clearing these vectors is not really needed
  index_v.clear();
  for(int i=0; i!= order+1; i++)
    {
      int v_size = m_Cov_xxt[i].m_size;// vertical size of the chunk under investigation
      for(int j=1; j!= 2*order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_Cov_xxt[j].m_size; // horizontal size of the chunk under investigation
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

		  double value 
		    =  m_Cov_xxt[i+j].find_monomial(monomial) 
		    +  m_xt[i].m_components[k] * m_x[j].m_components[l]
		    -  m_Cov_xxt[i].m_components[k] * m_x[j].m_components[l]
		    -  m_xt[i].m_components[k] * m_Cov_xx[j].m_components[l];


		  m_eigen_matrix(offset_xt[i]+k,size_xt + offset_xx[j]+l) = value;  
		  m_eigen_matrix(size_xt + offset_xx[j]+l,offset_xt[i]+k) = value;  

		  if(l!=h_size-1){
		    m_x[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_xt[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }


  index_v.resize(2*order+1); 
  index_h.clear();
  index_v.clear();
  for(int i=1; i!= 2*order+1; i++)
    {
      int v_size = m_Cov_xx[i].m_size;// vertical size of the chunk under investigation
      for(int j=1; j!= 2*order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_Cov_xx[j].m_size; // horizontal size of the chunk under investigation
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

		  double value 
		    =  m_Cov_xx[i+j].find_monomial(monomial) 
		    +  m_x[i].m_components[k] * m_x[j].m_components[l]
		    -  m_Cov_xx[i].m_components[k] * m_x[j].m_components[l]
		    -  m_x[i].m_components[k] * m_Cov_xx[j].m_components[l];


		  m_eigen_matrix(size_xt + offset_xx[i]+k,size_xt + offset_xx[j]+l) = value;  
		  m_eigen_matrix(size_xt + offset_xx[j]+l,size_xt + offset_xx[i]+k) = value;  

		  if(l!=h_size-1){
		    m_x[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_x[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }
  
  double const avg_weight = m_sum_weight/m_entries;
  m_eigen_matrix*=1.0/((T)(m_entries - 1)*avg_weight*avg_weight);
  //  return out;
}

/*
template <class T>
void full_correlation_tensor_serie<T>::cov_matrix_1dim(unsigned int order,unsigned int startorder)
{
  max_significant_order = 0;
  cov_vector(order);
  //std::cout << "order = " << order << std::endl;
  // std::vector<unsigned int> index_v(order+1,0);
  // std::vector<unsigned int> index_h(order+1,0);
  // std::vector<unsigned int> offset(order+1,0);
  // std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=startorder; i!=order+1; i++)
    {
      //      offset[i] = size;
      size += m_C[i].m_size;
    }
#ifdef USE_LAPACK
  m_matrix.resize(size);
#else
  m_eigen_matrix.resize(size,size);
#endif

  for(int i=startorder; i!= order; i++)
    {
      for(int j=i; j!= order; j++) // we fill a symmetric matrix
	{
	  double value =  m_C[i+j].m_components[0]- m_vector[i-startorder]*m_vector[j-startorder]; //vector is costructed by cov_vector: indexing skips one component
#ifdef USE_LAPACK
	  m_matrix( i , j ) = value; 
#else
	  m_eigen_matrix(i,j) = value;  
	  m_eigen_matrix(j,i) = value;  
#endif
	  //	  std::cout << "[" << i << "][" << j << "] = " << value << std::endl;
	}
    }
  m_eigen_matrix*=1.0/(T)m_entries;
  //  return out;
}
*/

template <class T>
void full_correlation_tensor_serie<T>::decompose(int order)
{
  //  normalize();
  //print("solver");
  //vector(m_max_F_order);
  if(m_dimensions == 1)
    {
      //matrix_1dim(m_max_h_order);
      matrix_1dim(order);
    }
  else
    {
      //matrix(m_max_h_order);
      //std::cout << "hi" << std::endl;
      matrix(order);
      //std::cout << "bye" << std::endl;
    }

  m_decomposed_matrix = m_eigen_matrix.ldlt();

}

template<class T>
tensor_serie_function<T> const full_correlation_tensor_serie<T>::solve(unsigned int const order)//decomposition & ext_decomposed_matrix)
{
  //  normalize();

  vector(order);
  unsigned int const n = m_vector.size();
  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(n);
  for(unsigned int i=0; i!=n ; i++)
    {v[i] = m_vector[i];}

  //decompose(order);
  
  std::cout << "matrix..." << std::endl;
  matrix(order);
  std::cout << "done "<< std::endl;
  Matrix Ginv = m_eigen_matrix.inverse();

  Eigen::Matrix<T, Eigen::Dynamic, 1> x = Ginv*v; //m_decomposed_matrix.solve(v);
  //m_vector.assign(&v.data()[0], &x.data()[n-1]);
 
  tensor_serie_function<T> func(m_dimensions, order);
  
  unsigned int offset = 0;
  for(int i=0; i!= order + 1; i++ )
    {
      unsigned int size = func.m_F[i].m_size;
      int k=0; 
      for(int j=offset; j!=offset+size; j++)
	{
	  func.m_F[i].m_components[k++] = x[j];
	}
      offset+=size;
    }
  unsigned int size_v = offset;
  std::vector<T>().swap(m_vector);

  double chi2 = -x.transpose()*v;

  func.m_chi2 = chi2;
  //=================
  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(2*order+1,0);
  std::vector<unsigned int> offsets(2*order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size_h = 0;
  for(int i=0; i!=2*order+1; i++)
    {
      offsets[i] = size_h;
      size_h += m_x[i].m_size;
    }

  std::cout << "cov matrix..." << std::endl;
  cov_matrix(order);
  std::cout << "done" << std::endl;
  double biaso = 0;
  for(int i=0; i!=size_v; i++){
    for(int j=0; j!=size_v; j++){
      std::cout << "Ginv("<<i<<","<<j<<") = "<<  Ginv(i,j) << " m(" <<i<<","<<j<<") = "<< m_eigen_matrix(i,j) << std::endl;

      biaso += Ginv(i,j)*m_eigen_matrix(i,j);
    }
  }
  std::cout << "order " << order << " biasO = " << biaso << std::endl; 

  double bias = 0;
  if(order>0) 
    {
      std::cout << "dGOverdgF matrix..."<< m_x.size() << std::endl;
      Matrix dGOverdgF(size_v, size_h - offsets[1]);
      dGOverdgF.setZero();
      
      for(int i=0; i!= order+1; i++)
	{
	  int v_size = m_x[i].m_size;// vertical size of the chunk under investigation
	  for(int j=0; j!= order+1; j++) // we fill a symmetric matrix

	    ///should be  for(int j=0; j!= order+1; j++) /// and access m_F[j].m_comp[k] when dGoverdg is 1

	    {
	      int h_size = m_x[j].m_size; // horizontal size of the chunk under investigation
	      //create a h_size*v_size matrix from m_tensor[i+j]
	      std::fill_n(index_v.begin(),i,0);
	      std::cout << "T0" << std::endl;
	      for(int k=0; k!= v_size; k++) 
		{
		  std::fill_n(index_h.begin(),j,0);
		  std::cout << "T00" << std::endl;
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
      		      std::cout << "T " << i+j <<  " " << monomial[0] << std::endl;
		      if(i+j==0) {continue;}
		      int horiz = offsets[j+i] - offsets[1] + m_x[i+j].index_from_monomial(monomial);
		      int verti = offsets[i]+k;
		      std::cout << "access " << verti << " , " << horiz <<  " " << i << "," << j << std::endl;
			
		      dGOverdgF(verti, horiz)+=func.m_F[j].m_components[l];
		      		      
		      std::cout << "OK" << std::endl;
		      if(l!=h_size-1){
			m_x[j].update_index(index_h);
			//std::cout << "h update" << std::endl;
		      }
		      std::cout << "T2" << std::endl;
		    }
		  if(k!=v_size-1){m_x[i].update_index(index_v);
		    //std::cout << "v update" << std::endl;
		  }
		  std::cout << "T3" << std::endl;
		}
	    }
	}
      std::cout << "done" << std::endl;
      Matrix d2Eoverdgdh = Ginv*dGOverdgF;
      Matrix d2Eoverdgdg = dGOverdgF.transpose()*d2Eoverdgdh;

      std::cout << "bias..." << std::endl;
      for(int i=0; i!=size_h - offsets[1]; i++)
	{
	  for(int j=0; j!=size_v; j++) 
	    {
	      bias+= -4*d2Eoverdgdh(j,i)*m_eigen_matrix(j,size_v+i);
	    }
	}
      for(int i=0; i!=size_h - offsets[1]; i++)
	{
	  for(int j=0; j!= size_h - offsets[1]; j++)
	    {
	      bias+= 3*d2Eoverdgdg(i,j)*m_eigen_matrix(size_v+i,size_v+j);
	    }
	}
      std::cout << "done" << std::endl;
    }
  

  std::cout << "check normalization!!" << std::endl;
 
  func.m_bias = bias+biaso;
  func.m_biaso = biaso;
  func.m_n_params = n;

  return func;

}

template <class T>
tensor_serie_function<T> const * full_correlation_tensor_serie<T>::principal_axis()
{
  //creating the matrices could be solved by calling matrix(1), but I want to leave m_eigen_matrix intact
  //must be called after normalization!
  if(m_max_order == 0)
    {
      std::cout << "principal_axis() : max_order ==0" << std::endl;
      return 0;
    }

  if(m_dimensions == 1)
    {
      tensor_serie_function<T> *pr_axis = new tensor_serie_function<T>(m_dimensions, 1);
      pr_axis->m_F[0].m_components[0] = -m_x[1].m_components[0];
      pr_axis->m_F[1].m_components[0] = 1;
      return pr_axis;
    }

  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(m_dimensions);
  for(int i=0; i!= m_dimensions; i++)
    {
      v[i] = m_x[1].m_components[i];
    }

  Matrix M; M.resize(m_dimensions, m_dimensions);
  int size = m_dimensions;
  int k = 0; //the components of the d*d second degree symmetric tensor are serialised, I can access them in a symple way
  for(int i=0; i!= m_dimensions; i++)
    {
      for(int j=i; j!= m_dimensions; j++)
	{
	  T value = m_x[2].m_components[k++];
	  M(i, j) = value;
	  M(j, i) = value;	  
	}
    }
  M -=  v*v.transpose();
  Eigen::SelfAdjointEigenSolver<Matrix > eigensolver(M);
  if(eigensolver.info()!= Eigen::Success)
    {
      std::cout << "principal_axis() : egensolver.info()!= Success" << std::endl;
      return 0;
    }
  //T max_eigenvalue = eigensolver.eigenvalues()[m_dimensions-1]; //M is positive semidefinite
  Eigen::Matrix<T, Eigen::Dynamic, 1> principal_axis_ev = eigensolver.eigenvectors().col(m_dimensions-1);

  tensor_serie_function<T> *pr_axis = new tensor_serie_function<T>(m_dimensions, 1);
  (*pr_axis).m_F[0].m_components[0] = -v.transpose()*principal_axis_ev;
  for(int i = 0; i!= m_dimensions; i++)
    {
      (*pr_axis).m_F[0].m_components[i] = principal_axis_ev[i];
    }

  if(pr_axis==0) {std::cout << "returning NULL pr_axis" << std::endl;}
  else   {std::cout << "returning valid pr_axis" << std::endl;}
  return pr_axis;
}



template <class T>
void full_correlation_tensor_serie<T>::print(std::string name)
{
  std::cout << "___" << std::endl;
  std::cout << "h Name = '"<< name << "' Max order = " << m_max_order << std::endl;
  for(int i=0 ; i!=m_max_order+1; i++)
    {
      m_xt[i].print();
    }  
  std::cout << "---" << std::endl; 
  for(int i=0 ; i!=2*m_max_order+1; i++)
    {
      m_x[i].print();
    }

  std::cout << m_eigen_matrix << std::endl;
  std::cout << "The determinant of A is " << m_eigen_matrix.determinant() << std::endl;
  std::cout << "The inverse of A is:\n" << m_eigen_matrix.inverse() << std::endl;
}

template <class T>
void full_correlation_tensor_serie<T>::release()
{
  std::vector<symmetric_tensor<T> >().swap(m_xt);
  std::vector<symmetric_tensor<T> >().swap(m_x);
  std::vector<symmetric_tensor<T> >().swap(m_Cov_xx);
  std::vector<symmetric_tensor<T> >().swap(m_Cov_xxt);
  std::vector<symmetric_tensor<T> >().swap(m_Cov_xtxt);
}


template <class T>
tensor_serie_function<T>::tensor_serie_function(tensor_serie_function<T> const & input)
{
  m_F = input.m_F;
  m_chi2 = input.m_chi2;
  m_bias = input.m_bias;
  m_biaso = input.m_biaso;
  m_n_params = input.m_n_params;
  m_F_derivate = input.m_F_derivate;
  m_Cov = input.m_Cov;
  m_max_order = input.m_max_order;
  m_dimensions = input.m_dimensions;
}

template <class T>
tensor_serie_function<T>::tensor_serie_function(unsigned int dimensions,unsigned int max_order)
{
  m_max_order = max_order;

  m_F.clear();
  m_F.resize(m_max_order+1);
  m_dimensions = dimensions;
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_F[i].init(m_dimensions, i);
    }
}


template <class T>
T const tensor_serie_function<T>::eval(tensor_serie<T> const & x ) const 
{
  T value = 0.0;
  for(int i=0; i!= m_max_order+1; i++)
    {
      if(i>0)
	{
	  value+=m_F[i].eval_diad(x.m_S[i].m_components);
	}
      else
	{
	  value+=m_F[i].m_components[0];
	}
    }
  return value;
}

template class tensor_serie_function<float>;
template class tensor_serie_function<double>;
template class tensor_serie_function<long double>;

template class full_correlation_tensor_serie<float>;
template class full_correlation_tensor_serie<double>;
template class full_correlation_tensor_serie<long double>;

