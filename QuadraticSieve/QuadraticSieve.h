#if !defined _QUADRATIC_SIEVE_H
#define _QUADRATIC_SIEVE_H

#include "bignum.h"

class QuadraticSieve{
public:
	BigNumber modPow(const BigNumber& a, const BigNumber& k, const BigNumber& n);
	BigNumber LegendreSymbol(BigNumber& a, BigNumber& p);
	QuadraticSieve(BigNumber& n);
	std::pair<BigNumber, BigNumber> doFactorization();
	vector<pair<BigNumber, vector<bool>>> sieving();
	BigNumber Tonelli_Shanks(BigNumber& a, BigNumber& p);
private:
	BigNumber N;
	std::pair<BigNumber, BigNumber> divisors;
	std::vector<BigNumber> Base;
	Ipp32u fbSize;
};

#endif