#include "bignum.h"
#include <string>
#include <map>

class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void checkPrime();
	void init();
	Ipp32u tDivBound = 1000;
	vector<Ipp32u> tDiv;
	std::map<BigNumber, Ipp32u> getFactor();
private:
	std::map<BigNumber, Ipp32u> factor;
	BigNumber number;
	bool isPrime = false;
	int nTraits = 6;
};