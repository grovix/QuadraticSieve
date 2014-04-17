#include "bignum.h"
#include <string>

class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void checkPrime();

	Ipp32u trialDivisionBound = 1000;
private:
	std::vector<BigNumber> factor;
	BigNumber number;
	bool isPrime = false;
	int nTraits = 6;
};