/*
 * variable_selection.cpp
 *
 *  Created on: 19 Aug 2020
 *      Author: kovesarki
 */



#include <Eigen/Dense>
#include "include/monopoly.h"
#include <random>
#include <iostream>
#include <vector>
#include <string>

#include "include/tensor_serie.hh"
#include "include/correlation_tensor_serie.hh"

int main(void) {
	std::vector<std::vector<double> > signal(10000, std::vector<double>(3,0.0) );
	std::vector<std::vector<double> > background(10000, std::vector<double>(3,0.0) );

	std::random_device rd;
	std::mt19937 gen(rd());

	std::normal_distribution<> sx(0,1.01);
	std::normal_distribution<> sy(1.0,0.4);
	std::normal_distribution<> sz(0.7,2.0);

	std::normal_distribution<> bx(0.5,0.9);
	std::normal_distribution<> by(1.4,0.7);
	std::normal_distribution<> bz(0.73,2.01);

	std::vector<std::string >  variable_names(3);
	variable_names[0] = "parameter x";
	variable_names[1] = "parameter y";
	variable_names[2] = "parameter z";

	for(int i=0; i!= signal.size(); i++ ){
		signal[i][0] = sx(gen);
		signal[i][1] = sy(gen);
		signal[i][2] = sz(gen);
	}

	for(int i=0; i!= background.size(); i++ ){
		background[i][0] = bx(gen);
		background[i][1] = by(gen);
		background[i][2] = bz(gen);
	}

	tensor_serie<double> m_x(3,2);
	correlation_tensor_serie<double> signal_fourier(3,1); //number of dimensions, number of degree (linear is enough for comparison)
	for(int i=0; i!= signal.size(); i++){
		m_x.create_diad(signal[i]);
		signal_fourier.fill(m_x,1.0,1.0); //training for target y=1.0,weights=1.0
	}
	signal_fourier.normalize();

	correlation_tensor_serie<double> background_fourier(3,1); //number of dimensions, number of degree (linear is enough for comparison)
	for(int i=0; i!= background.size(); i++) {
		m_x.create_diad(background[i]);
		background_fourier.fill(m_x,1.0,1.0);
	}
	background_fourier.normalize();

	signal_fourier.cov_matrix(1);
	background_fourier.cov_matrix(1);

	double a = signal_fourier.sigma_distance(background_fourier);
	std::cout << "significance is " << a << std::endl;

	std::vector<double> sigmas(variable_names.size(), 0.0);
	std::vector<int> variable_ordering(variable_names.size(), 0.0);

	signal_fourier.variable_selection(background_fourier, sigmas, variable_ordering);

	std::cout << "Variable order, from best to worst. " << std::endl
			  << "Significance is meant for the group of best variables" << std::endl
			  << "(the variable, and all the variables above)." << std::endl;
	for(int i=0; i!= variable_ordering.size(); i++)
	{
	  std::cout << sigmas[variable_ordering[i]] << " \t" << variable_names[variable_ordering[i]] << std::endl;
	}

	std::cout << "varnames["<< variable_ordering.size()<< "] = {" << std::endl;
	for(int i=variable_ordering.size()-1; i!= -1; i--)
	{
	  std::cout << "\"" << variable_names[variable_ordering[i]] << "\"";
	  if(i!=0) {std::cout << ","; }
	  else {std::cout << std::endl << "}";}
	  std::cout << std::endl;
	}
}

