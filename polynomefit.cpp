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

	tensor_series<double> polynom(1,10);

	for(int i=0; i!= data.size(); i++) {
		//cout<< i << data[i][0] << " " << data[i][1] << endl;

		polynom.fill_1dim(data[i][0], data[i][1], 1.0);
	}
	polynom.normalize();
	polynom.solve();

	for(double x=0; x< 1.0; x+=0.1) {
		double y = (-0.1*x*x*x + 0.7*x*x + 0.01*x + 0.5);
		cout << "p("<<x<<")=" << polynom.eval_1dim(x) << " vs " << y << endl;
	}

	//polynom.print("10 degree fit");



}
