#include "../include/fractional_fit.hh"
 //#include "Math/RootFinder.h"
#include "TMath.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/Minimizer.h"
//#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include <iostream>


template<class T>
fractional_fit<T>::fractional_fit(Matrix const & A, Matrix const & B,Matrix const & D,
				  //, Vector const & m, Vector const & n, Vector const & d)
				  std::vector<T> const &m, std::vector<T> const &n, std::vector<T> const &d)
{
  this->A = A;
  this->B = B;
  this->D = D;
  
  this->m.resize(m.size());
  this->n.resize(m.size());
  this->d.resize(m.size());

  for(int i=0; i!=m.size() ; i++){
    this->m[i] = m[i];    
    this->n[i] = n[i];
    this->d[i] = d[i];
  }
}



template<class T>
void fractional_fit<T>::analytic_solve()
{
  Vector diff = n-d;
  Vector t = (B+D)*((A-B).ldlt().solve(m-n));
  std::cout << " n-d = " << diff << std::endl;
  std::cout << " t = " << t << std::endl;

  std::cout << " a = " << t.norm()/diff.norm() << std::endl;

}


template<class T>
void fractional_fit<T>::solve()
{
  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  //   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("Minuit2", "Fumili");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLMultiFit", "");
//   ROOT::Math::Minimizer* min = 
//          ROOT::Math::Factory::CreateMinimizer("GSLSimAn", "");

  //min->SetMaxFunctionCalls(1000000);
  //min->SetMaxIterations(100000);
  //min->SetTolerance(0.001);
  //ROOT::Math::Functor f(&RosenBrock,2); 
  
  ROOT::Math::Functor1D func(*this);
  
  ROOT::Math::BrentMinimizer1D bm;
  bm.SetFunction(func, 0.0,1.0);
  bm.Minimize(1000,0.0001,0.0001);
  a_min = bm.XMinimum();
  chi2 = bm.FValMinimum();
}

template<class T>
double fractional_fit<T>::operator()(double a)
{
  double a2 = a*a;
  Sigma = a2*A + (1.0 - a2)*B + D;
  l = a*m + (1.0-a)*n - d;
  return l.transpose()*Sigma.ldlt().solve(l);
}

template class fractional_fit<float>;
template class fractional_fit<double>;
template class fractional_fit<long double>;
