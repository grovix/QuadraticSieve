#include "bignum.h"
#include <string>
#include <map>

class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void init();
	vector<BigNumber> tDiv;
	std::map<BigNumber, Ipp32u> getFactor();
	void rho_Pollard(BigNumber N);
	void QuadraticSieve(BigNumber N);
	Ipp32u tDivBound = 1000;
	Ipp32u pollardIter = 1000000;
	int nTraits = 100;
	void insertDivisor(BigNumber& a);
private:
	std::map<BigNumber, Ipp32u> factor;
	std::vector<BigNumber> intermNumbers;
	BigNumber number;
	bool isPrime = false;
	bool isFactoredByPollard = false;
};