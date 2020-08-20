//============================================================================
// Name        : polynomfit.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

//#include <stdio.h>
//#include <stdlib.h>
#include <Eigen/Dense>
#include "monopoly.h"
#include "include/tensor_serie.hh"
#include "full_correlation_tensor_serie.hh"
#include <random>

#include <iostream>

using namespace std;

int main(void) {
//	puts("Hello World!!!");
//	return EXIT_SUCCESS;

	cout << "Hello world" << endl;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic > m_eigen_matrix;
	m_eigen_matrix.resize(2,2);

	m_eigen_matrix << 1, 2,
					3,4;
	cout << m_eigen_matrix << endl;


	std::vector<std::vector<double> > data(10000, std::vector<double>(2,0.0) );
	std::random_device rd;
	std::mt19937 gen(rd());

	std::normal_distribution<> d(0,1);
	std::uniform_int_distribution<> distrib(1, 6);
	for(int i = 0; i!= data.size(); i++) {
		const double x =  distrib(gen) + d(gen);
		const double y = -0.1*x*x*x + 0.7*x*x + 0.01*x + 0.5 + d(gen);
		data[i][0] = x;
		data[i][1] = y;
	}

	full_correlation_tensor_serie<double> polynom(1,10);
	tensor_serie<double> m_x(1,10*4); //the internal tensors orders are 4 times larger

	for(int i=0; i!= data.size(); i++) {
		//cout<< i << data[i][0] << " " << data[i][1] << endl;
		m_x.create_diad_1dim(data[i][0]);
		polynom.fill_1dim(m_x, data[i][1], 1.0);
	}
	polynom.normalize();

	cout << "Solving for degree 4" << std::endl;
	tensor_serie_function<double> pol4 = polynom.solve(4);
	cout << "chi2\t expected chi2 \t the traditional bias \t full bias" << endl;
	cout << pol4.m_chi2 << "\t" << (pol4.m_chi2 + pol4.m_bias) << "\t" << pol4.m_biaso << "\t" << pol4.m_bias << endl;
	cout << "Solving for degree 10" << std::endl;
	tensor_serie_function<double> pol10 = polynom.solve(10);
	cout << "chi2\t expected chi2 \t the traditional bias \t full bias" << endl;
	cout << pol10.m_chi2 << "\t" << (pol10.m_chi2 + pol10.m_bias) << "\t" << pol10.m_biaso << "\t" << pol4.m_bias << endl;

	for(double x=0; x< 1.0; x+=0.1) {
		double y = (-0.1*x*x*x + 0.7*x*x + 0.01*x + 0.5);
		m_x.create_diad_1dim(x);
		cout << "p("<<x<<")=" << pol4.eval(m_x) << " vs " << y << endl;
	}

	//polynom.print("10 degree fit");



}
