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
#include "include/treepoly_smooth.hh"
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

	int n_levels = 5, n_input_vars =1;

	treepoly_smooth<long double> tree_poly(n_input_vars,n_levels);
	tree_poly.reserve(data.size());

	for(int i=0; i!= data.size(); i++) {
		//cout<< i << data[i][0] << " " << data[i][1] << endl;
		std::vector<double> xv(1);

		xv[0] = data[i][0];
		tree_poly.fill_1dim(data[i][0], data[i][1], 1.0);
	}
	tree_poly.train();

	for(double x=0; x< 1.0; x+=0.1) {
		double y = (-0.1*x*x*x + 0.7*x*x + 0.01*x + 0.5);
		std::vector<long double> xv(1);
		xv[0] = x;
		cout << "p("<<x<<")=" << tree_poly.eval(xv) << " vs " << y << endl;
	}

	//polynom.print("10 degree fit");



}
