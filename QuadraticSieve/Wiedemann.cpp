#include "wiedemann.h"
#include "bignum.h"

Wiedemann::Wiedemann(const SparseMatrix& a){
	B = a;
	N = B.getSize();
}

vector<bool> Wiedemann::getSolution(){
	vector<bool> x(getRandomVector());
	vector<bool> y(getRandomVector());
	vector<bool> z(B.Multiply(y));

	uInt L = 2 * N + delta;
	vector<bool> a(L + 2, false);
	a[0] = vecMultiply(x, z);

	vector<bool> Bt(B.Multiply(z));
	for (uInt i = 1; i < L+2; i++){
		a[i] = vecMultiply(x, Bt);
		Bt = B.Multiply(Bt);
	}

	vector<bool> T = Berlekamp_Massey(a);
	uInt l = T.size();
	vector<bool> F(l);
	std::copy(T.rbegin(), T.rend(), F.begin());
	T.clear();

	uInt p = 0;
	for (uInt i = 0; i < l; ++i){
		if (F[i]){
			p = i;
			break;
		}
	}
	vector<bool> w(y);
	Bt = B.Multiply(y);
	for (uInt i = p + 1; i < l; ++i){
		if (F[i])
			w = vecSum(w, Bt);
		Bt = B.Multiply(Bt);
	}

	x = w;
	w = B.Multiply(w);
	while (!B.isZero(w)){
		x = w;
		w = B.Multiply(w);
	}

	return w;
}

vector<bool> Wiedemann::getRandomVector(){

	BigNumber t;

	int size;
	// define Pseudo Random Generator (default settings)
	ippsPRNGGetSize(&size);
	IppsPRNGState* pPrng = (IppsPRNGState*)(new Ipp8u[size]);
	ippsPRNGInit(160, pPrng);
	// define 256-bits Big Numbers X and Y
	const int bnBitSize = 1;

	vector<bool> v(N);
	for (auto && i : v){
		ippsPRNGen_BN(BN(t), bnBitSize, pPrng);
		vector<Ipp32u> g;
		t.num2vec(g);
		i = (bool)g[0]; //Mother of Random !!!
	}
	return std::move(v);
}

bool Wiedemann::vecMultiply(vector<bool>& a, vector<bool>& b){

	uInt n = a.size();
	bool result = false;
	for (uInt i = 0; i < n; ++i){
		result^=(a[i] & b[i]);
	}
	return result;
}

vector<bool> Wiedemann::vecSum(vector<bool>& a, vector<bool>& b){

	vector<bool> res(a.size());
	for (uInt i = 0; i < a.size(); ++i){
		res[i] = a[i] ^ b[i];
	}
	return std::move(res);
}

vector<bool> Wiedemann::Berlekamp_Massey(vector<bool>& a){
	uInt M = a.size();
	vector<bool> b(M);
	vector<bool> C(M);
	vector<bool> T(M);
	vector<bool> F;

	b[0] = true;
	C[0] = true;
	uInt l = 0, m = -1;
	for (uInt n = 0; n < M; ++n){
		bool d = false;
		for (uInt i = 0; i <= l; ++i){
			d ^= C[i] & a[n - i];
		}
		if (d){
			std::copy(C.begin(), C.end(), T.begin());
	     	int Nm = n - m;
			for (uInt j = 0; j < M - Nm; ++j){
				C[Nm + j] = C[Nm + j] ^ b[j];
			}
			if (l <= n / 2){
				l = n + 1 - l;
				m = n;
				std::copy(T.begin(), T.end(), b.begin());
			}
		}
	}

	for (uInt i = 0; i < l + 1; ++i)
		F.push_back(C[i]);

	return std::move(F);
}