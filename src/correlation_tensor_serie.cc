#include "../include/correlation_tensor_serie.hh"
#include "../include/tensor_serie.hh"
#include <iostream>
//#include "../include/kissing_ellipsoid_solver.hh"
#include "../include/fractional_fit.hh"


template <class T>
correlation_tensor_serie<T>::~correlation_tensor_serie()
{
  //  std::cout << "correlation_tensor_serie::~correlation_tensor_serie() m_max_order = "<< m_max_order << std::endl;
}

template <class T>
correlation_tensor_serie<T>::correlation_tensor_serie(unsigned int dimensions,unsigned int max_order)
{
  m_max_order = max_order;

  m_S.clear();
  m_S.resize(m_max_order+1);
  m_dimensions = dimensions;
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_S[i].init(m_dimensions, i);
    }

  m_C.clear();
  m_C.resize(m_max_order*2+1);
  for(int i=0; i!=2*m_max_order+1; i++)
    {
      m_C[i].init(m_dimensions, i);
    }
  m_sum_weight = 0.0;
  m_entries=0;
}

template <class T>
void correlation_tensor_serie<T>::create_diad(std::vector<T> const & /*in_x*/ x)
{ 
#ifndef LOGDIAD
	    m_S[0].m_components[0] = 1;
#else
	    m_S[0].m_components[0] = 0;
#endif

  for(int i=0; i!=m_max_order+1; i++) 
    {
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_S[i].m_components = x;
#else
	    for(int j=0; j!=m_dimensions ; j++)
	      {
		m_S[i].m_components[j] = log(x[j]);
	      }
#endif
	  }
	  else
	    {
	      m_S[i].create_diad(m_S[1].m_components, m_S[i-1].m_components);	
	    }
	}
    }
}

template <class T>
void correlation_tensor_serie<T>::fill(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_S[i].fill_diad(x.m_S[i].m_components,target, weight );      
    }
  for(int i=0; i!=2*m_max_order+1; i++) 
    {
      m_C[i].fill_diad(x.m_S[i].m_components,target*target, weight*weight );      
    }
  m_sum_weight+=weight;
  m_entries++;
}

template <class T>
void correlation_tensor_serie<T>::create_diad_1dim(const T & x)
{
#ifndef LOGDIAD
	    m_S[0].m_components[0] = 1;
#else
	    m_S[0].m_components[0] = 0;
#endif
  for(int i=0; i!=2*m_max_order+1; i++) 
    {
      if(i>0)
	{
	  if(i==1) {
#ifndef LOGDIAD
	    m_S[i].m_components[0] = x;
#else
	    m_S[i].m_components[0] = log(x);
#endif
	  }
	  else
	    {
	      m_S[i].create_diad_1dim(m_S[1].m_components[0], m_S[i-1].m_components[0]);	
	    }
	}
    }
}
template <class T>
void correlation_tensor_serie<T>::fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_S[i].fill_diad_1dim(x.m_S[i].m_components[0],target, weight );
    }

  for(int i=0; i!=2*m_max_order+1; i++) 
    {
      m_C[i].fill_diad_1dim(x.m_S[i].m_components[0],target*target, weight*weight );
    }
  m_sum_weight+=weight;
  m_entries++;
}


template <class T>
void correlation_tensor_serie<T>::normalize()
{
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_S[i].normalize();      
    }
  for(int i=0; i!=2*m_max_order+1; i++)
    {
      m_C[i].normalize(m_sum_weight*m_sum_weight/(T)m_entries);      
    }
  //std::cout << "sum of weights = " << m_sum_weight << " entries = " << m_entries << std::endl;
}

template <class T>
void correlation_tensor_serie<T>::pretrain_add(correlation_tensor_serie<T> const & serie_b)
{
  for(int i=0; i!= m_max_order+1; i++)
    {
      m_S[i].add(serie_b.m_S[i]);      
    }
}

template <class T>
void correlation_tensor_serie<T>::vector(unsigned int order)
{
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      size += m_S[i].m_size;
    }
  m_vector.resize(size);
  //m_eigen_vector.resize(size);
  int k = 0;
  //std::cout << "size = " << size << std::endl;  
  for(int i=0; i!=order+1; i++)
    {
      for(int j=0; j!= m_S[i].m_size; j++)
	{
	  m_vector[k] = m_S[i].m_components[j];
	  k++;
	}
    }
  //std::cout << "end" << std::endl;
}

template <class T>
void correlation_tensor_serie<T>::cov_vector(unsigned int order, unsigned int startorder)
{
  unsigned int size = 0;
  for(int i=startorder; i!=order+1; i++) //this vector starts from i=1 
    {
      size += m_S[i].m_size;
    }
  m_vector.resize(size);
  //m_eigen_vector.resize(size);
  int k = 0;
  //std::cout << "size = " << size << std::endl;  
  for(int i=startorder; i!=order+1; i++)
    {
      for(int j=0; j!= m_S[i].m_size; j++)
	{
	  m_vector[k] = m_S[i].m_components[j];
	  k++;
	}
    }
  //std::cout << "end" << std::endl;
}


template <class T>
void correlation_tensor_serie<T>::matrix(unsigned int order)
{  
  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset(order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=0; i!=order+1; i++)
    {
      offset[i] = size;
      size += m_S[i].m_size;
    }
  
#ifdef USE_LAPACK
  m_matrix.resize(size);
  //  std::cout << "matrix size = " << m_matrix.get_size() << " " << size << std::endl;
#else
  m_eigen_matrix.resize(size, size);
#endif


  for(int i=0; i!= order+1; i++)
    {
      int v_size = m_S[i].m_size;// vertical size of the chunk under investigation
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_S[j].m_size; // horizontal size of the chunk under investigation
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
		  double value =  m_S[i+j].find_monomial(monomial);
		  
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
		    m_S[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_S[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }
  //  return out;
}

template <class T>
void correlation_tensor_serie<T>::matrix_1dim(unsigned int order)
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
      size += m_S[i].m_size;
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
	  double value =  m_S[i+j].m_components[0]; //every g has only one component
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
void correlation_tensor_serie<T>::cov_matrix(unsigned int order, unsigned int startorder)
{  
  //for density tensor series the first order is 1, because the normalization is not a free parameter, therefore it does not show up in the covariance matrix
  cov_vector(order);

  std::vector<unsigned int> index_v(order+1,0);
  std::vector<unsigned int> index_h(order+1,0);
  std::vector<unsigned int> offset(order+1,0);
  std::vector<unsigned int> monomial(m_dimensions,0);
  unsigned int size = 0;
  for(int i=startorder; i!=order+1; i++) 
    {
      offset[i] = size;
      size += m_C[i].m_size;
    }
  
#ifdef USE_LAPACK
  m_matrix.resize(size);
  //  std::cout << "matrix size = " << m_matrix.get_size() << " " << size << std::endl;
#else
  m_eigen_matrix.resize(size, size);
#endif


  for(int i=startorder; i!= order+1; i++)
    {
      int v_size = m_C[i].m_size;// vertical size of the chunk under investigation
      for(int j=i; j!= order+1; j++) // we fill a symmetric matrix
	{
	  int h_size = m_C[j].m_size; // horizontal size of the chunk under investigation
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
		  double value =  m_C[i+j].find_monomial(monomial) -  m_S[i].m_components[k]*m_S[j].m_components[l];
		  
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
		    m_C[j].update_index(index_h);
		    //std::cout << "h update" << std::endl;
		  }
		}
	      if(k!=v_size-1){m_C[i].update_index(index_v);
		//std::cout << "v update" << std::endl;
	      }
	    }
	}
    }
  m_eigen_matrix*=1.0/(T)m_entries;
  //  return out;
}


template <class T>
void correlation_tensor_serie<T>::cov_matrix_1dim(unsigned int order,unsigned int startorder)
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

template <class T>
//void correlation_tensor_serie<T>::set_values(TVectorT<double> const & vector)
void correlation_tensor_serie<T>::set_values()
{
  //m_S.clear();
  //m_S.resize(m_max_order+1);
  unsigned int offset = 0;
  //double const * p = vector.GetMatrixArray();
  //for(int i=0; i!=m_max_h_order+1; i++)
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_S[i].init(m_dimensions, i);
      unsigned int size = m_S[i].m_size;      
      //      m_tensors[i].m_components.assign(p+offset,p+offset+size);
      int k=0; 
      for(int j=offset; j!=offset+size; j++)
	{
	  //m_F[i].m_components[k++] = vector[j];//factorial[i];
	  m_S[i].m_components[k++] = m_vector[j];
	}
      offset+=size;      
    }
  std::vector<T>().swap(m_vector);
  for(int i=0; i!=m_max_order +1; i++)
    {
      //m_S[i].release_binomial();
    }
}

template <class T>
T correlation_tensor_serie<T>::eval(tensor_serie<T> const & x )
{
  T value = 0.0;
  for(int i=0; i!= m_max_order+1; i++)
    {
      if(i>0)
	{
	  value+=m_S[i].eval_diad(x.m_S[i].m_components);
	}
      else
	{
	  value+=m_S[i].m_components[0];
	}
    }
  return value;
}

template <class T>
T correlation_tensor_serie<T>::eval_1dim(tensor_serie<T> const & x)
{
  T value = 0.0;
  //  for(int i=0; i!= m_max_h_order+1; i++)
  for(int i=0; i!= m_max_order+1; i++)
    {
      if(i>0)
	{
	  value+=m_S[i].eval_diad_1dim(x.m_S[i].m_components[0]);
	}
      else
	{
	  value+=m_S[i].m_components[0]; //the 0th order has 1 element only
	}
    }
  return value;
}

template <class T>
void correlation_tensor_serie<T>::decompose(int order)
{
  normalize();
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

  //std::cout << "Solving with eigen" << std::endl;

  /*  int n = m_vector.size();//h_v_temp.GetNrows(); 
  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(n);
  for(int i=0; i!=n ; i++)
    {v[i] = m_vector[i];}
  Eigen::Matrix<T, Eigen::Dynamic, 1> x = 
    m_eigen_matrix.ldlt().solve(v);
  std::cout << "finished" << std::endl;

  for(int i=0; i!=n ; i++)
    {
      m_vector[i] = x[i];
    }
    set_values();*/

  m_decomposed_matrix = m_eigen_matrix.ldlt();
#endif
}

template<class T>
void correlation_tensor_serie<T>::solve(decomposition & ext_decomposed_matrix)
{
  normalize();
#ifdef USE_LAPACK
#else
  vector(m_max_order);
  int n = m_vector.size();
  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(n);
  for(int i=0; i!=n ; i++)
    {v[i] = m_vector[i];}
  Eigen::Matrix<T, Eigen::Dynamic, 1> x = ext_decomposed_matrix.solve(v);
  //m_vector.assign(&v.data()[0], &x.data()[n-1]);
 
  for(int i=0; i!=n ; i++)
     {
       m_vector[i] = x[i];
     }
  set_values();
#endif
}

/*template<class T>
void correlation_tensor_serie<T>::fractions(correlation_tensor_serie<T> const & A, 
					       correlation_tensor_serie<T> const & B, 
					       T & fraction, T & error)
{
  
  fractional_fit<T> fitter(A.m_eigen_matrix, B.m_eigen_matrix, m_eigen_matrix, A.m_vector, B.m_vector, m_vector  );

  fitter.solve();
  fraction = fitter.a_min;
  error = fitter.chi2; //will be the error later

  fitter.analytic_solve();

} */

template<class T>
void correlation_tensor_serie<T>::variable_selection(correlation_tensor_serie<T> const & ext_serie, std::vector<T> & sigmas, std::vector<int> & order)
{
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;

  int size =m_vector.size();
  Vector mu; mu.resize(m_vector.size());
  Vector nu; nu.resize(m_vector.size());

  for(int i=0; i!=m_vector.size() ; i++){
    mu[i] = m_vector[i];
    nu[i] = ext_serie.m_vector[i];
  }

  Matrix E = m_eigen_matrix + ext_serie.m_eigen_matrix;
  Vector l = mu - nu;
  T distance = l.norm();
  //Vector diff_norm = diff/distance;


  Vector c = 0*l;
  std::vector<int> inserted(size, 0);
  for(int i=0; i!= size; i++)
    {
      double s_curr = 0.0, s_max = 0.0;
      int include_me = -1;
      for(int j=0; j!=size; j++)
	{
	  if(inserted[j] == 1){continue;}
	  c[j] = l[j];
	  double s = c.transpose()*E*c;
	  double n = c.norm();
	  s = n*n/sqrt(s);
	  if(s>s_max)
	    {
	      s_max = s;
	      include_me = j;
	    }
	  c[j] = 0.0;
	}
      std::cout << "s = " << s_max << " i = " << i << " new = " << include_me << std::endl;                                                
      s_curr = s_max;
      c[include_me] = l[include_me];
      sigmas[include_me] = s_max;        
      order[i] = include_me;  
      inserted[include_me] = 1;
    }


  /* // upside down
  Vector c = l;
  double n = c.norm();
  double s_prev = c.transpose()*E*c;
  s_prev = n*n/sqrt(s_prev);


  std::vector< int > ignores(size, 0);
  std::cout << "size = " << size << std::endl;
  for(int i=0; i!= size; i++)
    {
      int dropme = -1;
      double smax = 0.0;
      
      for(int k=0; k!= size; k++)
	{
	  if(ignores[k]!=0) {continue;}
	  if(i!=size-1)
	    {
	      c[k] = 0.0;
	    }
	  double n = c.norm();
	  double s = c.transpose()*E*c;
	  s = n*n/sqrt(s);
	  //if(fabs(s-s_prev)<fabs(smax-s_prev))
	  if(s>smax)
	    {
	      smax = s;
	      dropme = k;
	    }
	  c[k] = l[k];
	}
    
      std::cout << "s = " << smax << " i = " << i << " dropme = " << dropme << std::endl;
      s_prev = smax;
      c[dropme] = 0.0;
      ignores[dropme] = 1;
      sigmas[dropme] = smax;
      order[i] = dropme;
      } */

  /*
  std::cout << "size = " << size << std::endl;
  std::vector< int > ignores(size, 0);
  for(int i=0; i!= size ; i++)
    {
      T est_s_var = 0;
      T norm = l.norm();
      T common_term = l.transpose()*E*l;
      
      T significance = TMath::Sqrt(common_term);
      common_term = 2.0/TMath::Power(norm, 5) *common_term;
      Vector ds = (2.0/TMath::Power(norm, 4) )*E*l;

      std::cout << "significance = " << norm*norm/significance << std::endl;
      int ignore_index = -1;
      double max = 0.0;
      for(int j=0; j!=size; j++)
	{
	  if(ignores[j]==0)
	    {
	      ds[j] = (ds[j] + common_term)*l[j];
	      //std::cout << "ds[" << j << "] = " << ds[j] << std::endl;
	      if(1.0/fabs(ds[j])>max )
		{
		  max = 1.0/fabs(ds[j]);
		  ignore_index = j;
		}
	    }	 
	}
      std::cout << "ignore j = " << ignore_index << " significance estimate = " << 1.0/max*TMath::Power(norm/significance, 5) << std::endl;
      l[ignore_index] = 0.0; 
      ignores[ignore_index]=1; 
      sigmas[ignore_index] = norm*norm/significance;
      order[i] = ignore_index;
    }

  */
  


}


template<class T>
T correlation_tensor_serie<T>::sigma_distance(correlation_tensor_serie<T> const & ext_serie)
{
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;

  Vector mu; mu.resize(m_vector.size());
  Vector nu; nu.resize(m_vector.size());

  for(int i=0; i!=m_vector.size() ; i++){
    mu[i] = m_vector[i];
    nu[i] = ext_serie.m_vector[i];
  }

  //Matrix Sigma = m_eigen_matrix + ext_serie.m_eigen_matrix;
  //Matrix Sigma = m_eigen_matrix.inverse() + ext_serie.m_eigen_matrix.inverse();
  
  /* //first attempt: convoluted-gauss sigma distance
  Vector difference  = mu-nu;
  T distance = difference.norm();
  std::cout << "distance = " << distance << std::endl;
  Vector normal = difference*1.0/distance;

  //T sigma = normal.transpose()*Sigma*normal;
  //std::cout << "sigma_1^-2 = " << sigma << std::endl;
  //sigma = 1.0/sigma;
  T sigma_1 = normal.transpose()*m_eigen_matrix.inverse()*normal;
  T sigma_2 = normal.transpose()*ext_serie.m_eigen_matrix.inverse()*normal;

  T sigma = 1.0/sigma_1 + 1.0/sigma_2;
  std::cout << "sigma_1^-2 = " << sigma_1 << std::endl;
  std::cout << "sigma_2^-2 = " << sigma_2 << std::endl;
  std::cout << "sigma^2 = " << sigma << std::endl;
  //sigma = 1.0/sigma;
  return distance/TMath::Sqrt(sigma);
  */
  
  //kissing elipsoids distance

  //std::cout << "S_1:" <<  m_eigen_matrix << std::endl;;
  //std::cout << "S_2:" <<  ext_serie.m_eigen_matrix << std::endl;;
  
  Matrix E = m_eigen_matrix + ext_serie.m_eigen_matrix;
  Vector diff = mu - nu;
  T distance = diff.norm();
  Vector diff_norm = diff/distance;

  
  T sigma = diff_norm.transpose()*E*diff_norm;
  std::cout << "sigma = " << sigma << std::endl;
  
  return distance/sqrt(sigma);

  //kissing_ellipsoid_solver<T> helper(m_eigen_matrix, ext_serie.m_eigen_matrix, mu, nu);
  //return helper.solve();

  // Vector x = (m_eigen_matrix.inverse() - ext_serie.m_eigen_matrix.inverse()).inverse() * (m_eigen_matrix.inverse()*mu - ext_serie.m_eigen_matrix.inverse()*nu);
  // //Vector discrim = 
  // std::cout << "x  = " << x << std::endl;

  // T sigma_1 = (x - mu).transpose()*m_eigen_matrix.inverse()*(x - mu);
  // T sigma_2 = (x - nu).transpose()*ext_serie.m_eigen_matrix.inverse()*(x - nu);
  // std::cout << "sigma_1 = " << sigma_1 << " should be equal to sigma_2 = " << sigma_2 << std::endl;
  // return TMath::Sqrt(sigma_1);

  //this part calculates the overlap of the two: the p-value of occurance
  // Matrix Sigma = m_eigen_matrix.inverse() + ext_serie.m_eigen_matrix.inverse();
  // Eigen::Matrix<T, Eigen::Dynamic, 1> kappa = m_eigen_matrix.inverse() * mu +  ext_serie.m_eigen_matrix.inverse() * nu;
  // T exponent = mu.transpose() * m_eigen_matrix * mu +
  //   nu.transpose() * ext_serie.m_eigen_matrix * nu  -
  //   kappa.transpose() * Sigma.inverse() * kappa;
  
  // T probability = 2.0*TMath::Pi()*m_eigen_matrix.determinant()*ext_serie.m_eigen_matrix.determinant()*Sigma.determinant();
  // probability = TMath::Power(probability, - (T)(m_max_order+1)/2.0  ) *exp(-0.5*exponent);
  // return probability;
  
}



template <class T>
void correlation_tensor_serie<T>::print(std::string name)
{
  std::cout << "___" << std::endl;
  std::cout << "h Name = '"<< name << "' Max order = " << m_max_order << std::endl;
  for(int i=0 ; i!=m_max_order+1; i++)
    {
      m_S[i].print();
    }  
  std::cout << "---" << std::endl; 
  for(int i=0 ; i!=2*m_max_order+1; i++)
    {
      m_C[i].print();
    }

  std::cout << m_eigen_matrix << std::endl;
  std::cout << "The determinant of A is " << m_eigen_matrix.determinant() << std::endl;
  std::cout << "The inverse of A is:\n" << m_eigen_matrix.inverse() << std::endl;
}

template <class T>
void correlation_tensor_serie<T>::release()
{
  std::vector<symmetric_tensor<T> >().swap(m_S);
  std::vector<symmetric_tensor<T> >().swap(m_C);
}

template class correlation_tensor_serie<float>;
template class correlation_tensor_serie<double>;
template class correlation_tensor_serie<long double>;
