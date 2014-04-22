#if !defined _QUADRATIC_SIEVE_H
#define _QUADRATIC_SIEVE_H

#include "bignum.h"

class QuadraticSieve{
public:
	BigNumber modPow(const BigNumber& a, const BigNumber& k, const BigNumber& n);
	BigNumber LegendreSymbol(BigNumber& a, BigNumber& p);
	QuadraticSieve(BigNumber& n);
	std::pair<BigNumber, BigNumber> doFactorization();
private:
	BigNumber N;
	std::pair<BigNumber, BigNumber> divisors;
	std::vector<BigNumber> Base;
};

#endif