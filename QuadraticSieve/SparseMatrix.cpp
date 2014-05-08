#include "sparse_matrix.h"

SparseMatrix::SparseMatrix(){
}

SparseMatrix::SparseMatrix(const SparseMatrix& m){
	R = m.R;
	n = m.n;
}

SparseMatrix::SparseMatrix(std::vector<uInt> m, uInt size){
	R = m;
	n = size;
}

vector<uInt> SparseMatrix::getR(){
	return R;
}

uInt SparseMatrix::getSize(){
	return n;
}

vector<bool> SparseMatrix::Multiply(vector<bool> v){
	uInt ind = 0;
	vector<bool> res(n,false);
	if (isZero(v))
		return std::move(res);
	for (uInt i = 0; i < n; ++i){
		uInt len = R[ind];
		++ind;
		bool temp = false;
		for (uInt j = 0; j < len; ++j){
			temp = v[R[ind]] ^ temp;
			++ind;
		}
		res[i] = temp;
	}
	return std::move(res);
}

bool SparseMatrix::isZero(const vector<bool> v){
	bool res = true;
	for (const auto& i : v)
		if (i != 0){
			res = false;
			break;
		}
	return res;
}