#include "Factorization.h"

Factorization::Factorization(const char*s){
	number = BigNumber(s);
}

Factorization::Factorization(Ipp32u value){
	number = BigNumber(value);
}

void Factorization::checkPrime(){
	Factorization::isPrime = number.isPrime(Factorization::nTraits);
}