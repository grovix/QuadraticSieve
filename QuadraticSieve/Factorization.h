#if !defined _FACTORIZATION_H
#define _FACTORIZATION_H

#include "bignum.h"
#include "QuadraticSieve.h"
#include <string>
#include <map>

class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void init();
	vector<BigNumber> tDiv;
	std::map<BigNumber, Ipp32u> getFactor();
	void rho_Pollard(BigNumber& N);
	void CallQuadraticSieve(BigNumber& N);
	void insertDivisor(const BigNumber& a);
	void perfectPowerTest(BigNumber& a);

	Ipp32u tDivBound = 50000;
	Ipp32u pollardIter = 10;
	int nTraits = 10;
private:
	std::map<BigNumber, Ipp32u> factor;
	BigNumber number;
	bool isPrime = false;
	bool isFactoredByPollard = false;
};

#endif