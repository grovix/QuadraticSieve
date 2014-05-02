#if !defined _QUADRATIC_SIEVE_H
#define _QUADRATIC_SIEVE_H

#include "bignum.h"
#include <fstream>;

class QuadraticSieve{
public:
	BigNumber modPow(const BigNumber& a, const BigNumber& k, const BigNumber& n);
	BigNumber LegendreSymbol(BigNumber& a, BigNumber& p);
	QuadraticSieve(BigNumber& n);
	std::pair<BigNumber, BigNumber> doFactorization();
	vector<pair<BigNumber, vector<Ipp32u>>> sieving();
	BigNumber Tonelli_Shanks(BigNumber& a, BigNumber& p);
	BigNumber Q(BigNumber& x);
	int nTrials = 10;
private:
	BigNumber N;
	std::pair<BigNumber, BigNumber> divisors;
	std::vector<BigNumber> Base;
	Ipp32u fbSize;
	BigNumber A, B, C;
};

#endif