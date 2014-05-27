#if !defined _QUADRATIC_SIEVE_H
#define _QUADRATIC_SIEVE_H

#include "bignum.h"
#include <fstream>;
#include <omp.h>
#include <tbb\blocked_range.h>
#include <tbb\parallel_for.h>
typedef vector<pair<BigNumber, vector<Ipp32u>>> vPair;

class QuadraticSieve{
public:
	BigNumber modPow(const BigNumber& a, const BigNumber& k, const BigNumber& n);
	BigNumber LegendreSymbol(const BigNumber& a, const BigNumber& p);
	QuadraticSieve(BigNumber& n);
	std::pair<BigNumber, BigNumber> doFactorization();
	vPair sieving();
	BigNumber Tonelli_Shanks(const BigNumber& a,const BigNumber& p);
	BigNumber Q(const BigNumber& x);
	int nTrials = 10;
	vector<unsigned int> getSparseMatrix(vPair& v);
private:
	BigNumber N;
	std::vector<BigNumber> Base;
	Ipp32u fbSize;
	BigNumber A, B, C;
};

#endif