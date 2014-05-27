#if !defined _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H
#include <iostream>
#include <vector>
#include <iterator>
#include <omp.h>
#include <tbb\blocked_range.h>
#include <tbb\parallel_for.h>
typedef unsigned int uInt;
using namespace std;

class SparseMatrix{
public:
	SparseMatrix();
	SparseMatrix(const SparseMatrix& m);
	SparseMatrix(std::vector<uInt>& m, uInt size);
	vector<bool> Multiply(vector<bool>& v);
	vector<uInt> getR();
	uInt getSize();
	bool isZero(const vector<bool>& v);
private:
	uInt n;
	std::vector<uInt> R;
};

#endif