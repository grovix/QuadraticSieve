#include "Factorization.h"

Factorization::Factorization(std::string& value){
	number = BigNumber(value.c_str());
}

Factorization::Factorization(Ipp32u value){
	number = BigNumber(value);
}

void Factorization::checkPrime(){

}