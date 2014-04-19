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
	cout << "Trial division completed" << endl;

	this->checkPrime();

	if (!isPrime)
		rho_Pollard();
}

void Factorization::rho_Pollard(){
	int size;
	int numSize = this->number.BitSize();
	int counter = this->pollardIter;
	ippsPRNGGetSize(&size);
	IppsPRNGState* pPrng = (IppsPRNGState*)(new Ipp8u[size]);
	ippsPRNGInit(160, pPrng);

	BigNumber x;
	ippsPRNGen_BN(BN(x), numSize, pPrng);
	while (x > number)
		ippsPRNGen_BN(BN(x), numSize, pPrng);
	BigNumber z(x);
	BigNumber p;
	bool flag = false;
	while (!flag && counter > 0){
		x = (x*x + BigNumber::One()) % number;
		z = (z*z + BigNumber::One()) % number;
		z = (z*z + BigNumber::One()) % number;
		p = (z - x).b_gcd(number);
		if (p > BigNumber::One())
			flag = true;
		--counter;
	}
	if (flag){
		if (p.isPrime(nTraits)){
			factor.insert(std::pair<BigNumber, Ipp32u>(p, 1));
			cout << "Rho-Pollard method has found prime factor of number!" << endl;
		}
		else{
			intermNumbers.push_back(p);
			cout << "Rho-Pollard method has found non-prime factor of number!" << endl;
		}
		number /= p;
	}
	else
		cout << "Rho-Pollard method has not found any factors of number!" << endl;
}

void Factorization::router(){

}
