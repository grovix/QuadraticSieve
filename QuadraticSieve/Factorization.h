#include "bignum.h"
#include <string>
class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void checkPrime();
private:
	std::vector<BigNumber> factor;
	BigNumber number;
	bool isPrime = false;
	int nTraits = 6;
	IppBitSupplier rndFunc;
	void* pRndParam;
};