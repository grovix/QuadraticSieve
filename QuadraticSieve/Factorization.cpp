#include "Factorization.h"

Factorization::Factorization(const char*s){
	number = BigNumber(s);
}

Factorization::Factorization(Ipp32u value){
	number = BigNumber(value);
}

std::map<BigNumber, Ipp32u> Factorization::getFactor(){
	init();
	return factor;
}

void Factorization::init(){
	if (number.isPrime(nTraits)){
		cout << "Input number is prime!" << endl;
		insertDivisor(number);
		return;
	}
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
			tDiv.push_back(BigNumber((Ipp32u)(i - a_begin)));
	}

	auto& f_end = factor.end();
	for (auto&& i : tDiv){
		while (number % i == BigNumber::Zero()){
			insertDivisor(BigNumber(i));
		}
	}
	cout << "Trial division completed" << endl;
	if (!number.isPrime(nTraits))
		rho_Pollard(number);
	else{
		if (number != BigNumber::One())
			insertDivisor(number);
	}
	return;
}

void Factorization::rho_Pollard(BigNumber& N){
	int size;
	int numSize = N.BitSize();
	auto counter = this->pollardIter;
	ippsPRNGGetSize(&size);
	IppsPRNGState* pPrng = (IppsPRNGState*)(new Ipp8u[size]);
	ippsPRNGInit(160, pPrng);

	BigNumber x;
	ippsPRNGen_BN(BN(x), numSize, pPrng);
	while (x > N)
		ippsPRNGen_BN(BN(x), numSize, pPrng);
	BigNumber z(x);
	BigNumber p;
	bool flag = false;
	while (!flag && counter > 0){
		x *= x;
		x += BigNumber::One();
		x %= N;
		z *= z;
		z += BigNumber::One();
		z %= N;
		z *= z;
		z += BigNumber::One();
		z %= N;
		p = (z - x).b_gcd(N);
		if (p > BigNumber::One())
			flag = true;
		--counter;
	}
	if (flag){
		if (N % p != BigNumber::Zero())
			throw("Error! Wrong rho-pollard result");
		N /= p;
		if (p.isPrime(nTraits)){
			insertDivisor(p);
			cout << "Rho-Pollard method has found prime factor of " << N << endl;
		}
		else{
			cout << "Rho-Pollard method has found non-prime factor of" << N << endl;
			rho_Pollard(p);
		}

		if (N.isPrime(nTraits))
			insertDivisor(N);
		else
			rho_Pollard(N);
	}
	else{
		cout << "Rho-Pollard method has not found any factors of " << N << endl;
		CallQuadraticSieve(N);
	}
}

void Factorization::insertDivisor(BigNumber& a){
	auto& f = factor.find(a);
	if (f != factor.end())
		f->second = f->second + 1;
	else
		factor.insert(std::pair<BigNumber, Ipp32u>(a,1));
	number /= a;
}

void Factorization::CallQuadraticSieve(BigNumber& N){
	//QuadraticSieve q(N);
	//std::pair<BigNumber, BigNumber> res = q.doFactorization();
	return;
}

void Factorization::perfectPowerTest(BigNumber& a){

}