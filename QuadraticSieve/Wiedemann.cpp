#include "wiedemann.h"

Wiedemann::Wiedemann(const SparseMatrix& a){
	B = a;
	N = B.getSize();
}

vector<bool> Wiedemann::getSolution(){
	vector<bool> w(N,false);

	vector<bool> x(getRandomVector());
	vector<bool> y(getRandomVector());
	vector<bool> z(B.Multiply(y));

	uInt L = 2 * N + delta;
	vector<bool> a(L + 2, false);
	vector<vector<bool>> pp;


	return w;
}

vector<bool> Wiedemann::getRandomVector(){
	vector<bool> v(N);
	for (auto && i : v){
		i = rand() % 2;
	}
	return std::move(v);
}