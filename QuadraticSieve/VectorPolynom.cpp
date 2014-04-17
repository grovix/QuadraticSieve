#include "vector_polynom.h"

VectorPolynom::VectorPolynom(int degree){
	vector<bool> poly(degree + 1);
	std::move(poly.begin(), poly.end(), polynom.begin()); //TODO: read about initialization of vector and constructors
	p_degree = degree;
}

VectorPolynom::VectorPolynom(vector<bool> poly){
	std::copy(poly.begin(), poly.end(), polynom.begin());
	p_degree = poly.size() -1;
}

