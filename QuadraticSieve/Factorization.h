#include "bignum.h"
#include <string>
class Factorization{
public:
	Factorization(Ipp32u value);
	Factorization(std::string& value);
	void checkPrime();
private:
	std::vector<BigNumber> factor;
	BigNumber number;
	bool isFactorized = false;
};