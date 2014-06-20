#if !defined _QUADRATIC_SIEVE_H
#define _QUADRATIC_SIEVE_H

#include "bignum.h"
#include <fstream>;
#include <omp.h>
#include <ipp.h>
//#include <tbb\blocked_range.h>
//#include <tbb\parallel_for.h>
typedef vector<pair<BigNumber, vector<Ipp32u>>> vPair;

class QuadraticSieve{
public:
	BigNumber modPow(const BigNumber& a, const BigNumber& k, const BigNumber& n);
	BigNumber LegendreSymbol(const BigNumber& a, const BigNumber& p);
	QuadraticSieve(BigNumber& n);
	std::pair<BigNumber, BigNumber> doFactorization();
	vPair sieving();
	BigNumber Tonelli_Shanks(const BigNumber& a, const BigNumber& p, BigNumber&m, BigNumber&b,
		BigNumber& q_TS, BigNumber& n_TS, BigNumber& z_TS, BigNumber&  two_TS, BigNumber& t_TS, BigNumber& y_TS, BigNumber& r_TS, BigNumber& x_TS, 
		Ipp32u e_TS, int thread_id);
	BigNumber Q(const BigNumber& x, BigNumber& A, BigNumber& B, BigNumber& C);
	int nTrials = 10;
	vector<unsigned int> getSparseMatrix(vPair& v);
	int numThreads = 4;
	IppsPRNGState** pPrng_TS = nullptr;
	static ofstream out;
	void init();
	int numSize_TS = 10;  //What size "n" should be?
private:
	BigNumber N;
	std::vector<BigNumber> Base;
	Ipp32u fbSize;
};

#endif