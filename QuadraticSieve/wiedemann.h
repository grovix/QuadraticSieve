#if !defined _WIEDEMANN_H
#define _WIEDEMANN_H
#include "sparse_matrix.h"

class Wiedemann{
public:
	Wiedemann(const SparseMatrix& a);
	vector<bool> getSolution();
private:
	SparseMatrix A;
	uInt N;
};

#endif