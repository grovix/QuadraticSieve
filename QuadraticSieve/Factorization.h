#include "bignum.h"
#include <string>
#include <map>

class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(const char*s);
	void checkPrime();
	void init();
	vector<Ipp32u> tDiv;
	std::map<BigNumber, Ipp32u> getFactor();
	void rho_Pollard();
	void router();
	Ipp32u tDivBound = 1000;
	Ipp32u pollardIter = 1000;
private:
	std::map<BigNumber, Ipp32u> factor;
	std::vector<BigNumber> intermNumbers;
	BigNumber number;
	bool isPrime = false;
	int nTraits = 6;
};