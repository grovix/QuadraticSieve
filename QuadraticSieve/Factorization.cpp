#include "Factorization.h"

Factorization::Factorization(const char*s){
	number = BigNumber(s);
}

Factorization::Factorization(Ipp32u value){
	number = BigNumber(value);
}

void Factorization::checkPrime(){
	isPrime = number.isPrime(nTraits,rndFunc,pRndParam);
	cout << "is prime" << isPrime << endl;
}