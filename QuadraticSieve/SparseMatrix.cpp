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

}

bool SparseMatrix::isZero(vector<bool> v){
	
}