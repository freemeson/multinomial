# multinomial
A library to store symmetric tensors and covariants of multivariate polynomial, with fast evaluation and random access.
Added functionality to make polynomial regression of multivariate data. The multivariade diad tensors
of data vectors are created when the fill(data-vector, target) function is called, and aggregated into the symmetric tensor
elements. 
The symmetic tensors can be used independently of this use case, and have very fast random access, with very small footprint as only the necessary components are stored. The multiplicity of an element is not calculated, as this component falls out in the regression equations, but it can be easily implemented.
Before usage, one must create a singleton of the binomial_singleton class, that stores the components of a binomial function.
The limit for a 64bit implementation from numerical errors appear around 20 degrees for double and 25 degrees for long double 
precison. Regression limits are dominated by the size of the matrix, due to inversion errors around sizes 2000x2000. 

The implementation uses the Eigen library, from http://eigen.tuxfamily.org, but it was tested with LAPACK as well.
For crosschecking the accuracy of the binomial calculation it also uses the binomial function from ROOT, https://root.cern.ch/ . 

Some descripton and example:

Polynomial expansion of the binary classification function
https://arxiv.org/abs/1203.5647

Multivariate regression and fit function uncertainty
https://arxiv.org/abs/1310.1022

For any questions, mail me at freemeson@googlemail.com

Peter Kovesarki
