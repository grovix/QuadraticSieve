#include "wiedemann.h"

Wiedemann::Wiedemann(const SparseMatrix& a){
	A = a;
	N = A.getSize();
}

vector<bool> Wiedemann::getSolution(){
	vector<bool> result;

	vector<bool> x(N);
	vector<bool> y(N);
	vector<bool> z(N);


	return result;
}