#if !defined _WIEDEMANN_H
#define _WIEDEMANN_H
#include "sparse_matrix.h"
#include <fstream>

class Wiedemann{
public:
	Wiedemann(const SparseMatrix& a);
	vector<bool> getSolution();
	vector<bool> getRandomVector();
	uInt delta = 5;
	bool vecMultiply(vector<bool>& a, vector<bool>& b);
	vector<bool> vecSum(vector<bool>& a, vector<bool>& b);
	vector<bool> Berlekamp_Massey(vector<bool>& a);
private:
	uInt N;
	SparseMatrix B;
	const int bnBitSize = 10;
	ofstream out;
};

#endif