#include "Factorization.h"

Factorization::Factorization(const char*s){
	number = BigNumber(s);
}

Factorization::Factorization(Ipp32u value){
	number = BigNumber(value);
}

void Factorization::checkPrime(){
	Factorization::isPrime = number.isPrime(Factorization::nTraits);
	if (isPrime && number != BigNumber::One()){
		factor.insert(std::pair<BigNumber, Ipp32u>(number, 1));
		cout << "Factorization completed!" << endl;
	}
}

std::map<BigNumber, Ipp32u> Factorization::getFactor(){
	return factor;
}

void Factorization::init(){
	//sieve of Eratosphenes for trial division
	vector<bool> arr(tDivBound + 1,true);
	auto& a_begin = arr.begin();
	auto& a_end = arr.end();
	Ipp32u ind, p;
	arr[0] = false;
	arr[1] = false;
	auto& bound = a_begin + sqrt(a_end - a_begin);
	for (auto& i = arr.begin() + 2; i <= bound; ++i){
		if (*i == true){
			p = i - a_begin;
			ind = p * p;
			while (ind < (a_end - a_begin)){
				arr[ind] = false;
				ind += p;
			}
		}
	}
	for (auto& i = arr.begin() + 2; i != a_end; ++i){
		if (*i)
			tDiv.push_back(i - a_begin);
	}

	auto& f_end = factor.end();
	for (auto&& i : tDiv){
		while (number % i == BigNumber::Zero()){
			auto& find = factor.find(i);
			if (find != f_end){
				find->second += 1;
			}
			else
				factor.insert(std::pair<BigNumber, Ipp32u>(i, 1));
			number /= i;
		}
	}

	this->checkPrime();
}

