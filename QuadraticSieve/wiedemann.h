#if !defined _WIEDEMANN_H
#define _WIEDEMANN_H
#include "sparse_matrix.h"

class Wiedemann{
public:
	Wiedemann(const SparseMatrix& a);
	vector<bool> getSolution();
	vector<bool> getRandomVector();
	uInt delta = 20;
private:
	SparseMatrix B;
	uInt N;
};

#endif