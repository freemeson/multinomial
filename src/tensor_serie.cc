#include "../include/tensor_serie.hh"
#include <iostream>

template <class T>
tensor_serie<T>::~tensor_serie()
{
  //  std::cout << "tensor_serie::~tensor_serie() m_max_order = "<< m_max_order << std::endl;
}

template <class T>
tensor_serie<T>::tensor_serie(unsigned int dimensions,unsigned int max_order)
{
  m_max_order = max_order;

  m_S.clear();
  m_S.resize(m_max_order+1);
  m_dimensions = dimensions;
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_S[i].init(m_dimensions, i);
    }
}

template <class T>
void tensor_serie<T>::create_diad(std::vector<T> const & /*in_x*/ x)
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
void tensor_serie<T>::fill(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_S[i].fill_diad(x.m_S[i].m_components,target, weight );      
    }
}

template <class T>
void tensor_serie<T>::create_diad_1dim(const T & x)
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
void tensor_serie<T>::fill_1dim(tensor_serie<T> const & x, T const & target, T const & weight)
{
  for(int i=0; i!=m_max_order+1; i++) 
    {
      m_S[i].fill_diad_1dim(x.m_S[i].m_components[0],target, weight );
    }
}


template <class T>
void tensor_serie<T>::normalize()
{
  for(int i=0; i!=m_max_order+1; i++)
    {
      m_S[i].normalize();      
    }
}

template <class T>
void tensor_serie<T>::pretrain_add(tensor_serie<T> const & serie_b)
{
  for(int i=0; i!= m_max_order+1; i++)
    {
      m_S[i].add(serie_b.m_S[i]);      
    }
}

template <class T>
void tensor_serie<T>::vector(unsigned int order)
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
void tensor_serie<T>::matrix(unsigned int order)
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
void tensor_serie<T>::matrix_1dim(unsigned int order)
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
//void tensor_serie<T>::set_values(TVectorT<double> const & vector)
void tensor_serie<T>::set_values()
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
T tensor_serie<T>::eval(tensor_serie<T> const & x ) const
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
T tensor_serie<T>::eval_1dim(tensor_serie<T> const & x) const
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
void tensor_serie<T>::decompose(int order)
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
  //m_decomposed_matrix = m_eigen_matrix.fullPivLu();
#endif
}

template<class T>
void tensor_serie<T>::solve(decomposition & ext_decomposed_matrix)
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


template <class T>
void tensor_serie<T>::print(std::string name) const 
{
  std::cout << "___" << std::endl;
  std::cout << "h Name = '"<< name << "' Max order = " << m_max_order << std::endl;
  for(int i=0 ; i!=m_max_order+1; i++)
    {
      m_S[i].print();
    }
  std::cout << "---" << std::endl; 
}

template <class T>
void tensor_serie<T>::release()
{
  m_S.clear();
}

template <class T>
tensor_serie<T> const * tensor_serie<T>::principal_axis()
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
      tensor_serie<T> *pr_axis = new tensor_serie<T>(m_dimensions, 1);
      pr_axis->m_S[0].m_components[0] = -m_S[1].m_components[0];
      pr_axis->m_S[1].m_components[0] = 1;
      return pr_axis;
    }

  Eigen::Matrix<T, Eigen::Dynamic, 1> v; v.resize(m_dimensions);
  for(int i=0; i!= m_dimensions; i++)
    {
      v[i] = m_S[1].m_components[i];
    }

  Matrix M; M.resize(m_dimensions, m_dimensions);
  int size = m_dimensions;
  int k = 0; //the components of the d*d second degree symmetric tensor are serialised, I can access them in a symple way
  for(int i=0; i!= m_dimensions; i++)
    {
      for(int j=i; j!= m_dimensions; j++)
	{
	  T value = m_S[2].m_components[k++];
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

  tensor_serie<T> *pr_axis = new tensor_serie<T>(m_dimensions, 1);
  (*pr_axis).m_vector.resize(m_dimensions+1);
  (*pr_axis).m_vector[0] = -v.transpose()*principal_axis_ev;
  for(int i = 0; i!= m_dimensions; i++)
    {
      (*pr_axis).m_vector[i+1] = principal_axis_ev[i];
    }
  (*pr_axis).set_values();
  if(pr_axis==0) {std::cout << "returning NULL pr_axis" << std::endl;}
  else   {std::cout << "returning valid pr_axis" << std::endl;}
  return pr_axis;
}


template class tensor_serie<float>;
template class tensor_serie<double>;
template class tensor_serie<long double>;
