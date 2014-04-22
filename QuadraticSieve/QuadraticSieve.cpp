#include "QuadraticSieve.h"
#include <memory>

QuadraticSieve::QuadraticSieve(BigNumber& n){
	N = n;
}

BigNumber QuadraticSieve::modPow(const BigNumber&a, const BigNumber& k, const BigNumber& n){
	int deg = k.BitSize();
	BigNumber b(BigNumber::One());
	if (k == BigNumber::Zero())
		return b;

	string s;
	k.num2hex(s);
	Ipp32u size = 0;
	Ipp32u pl = 3;
	while (s[pl] == '0')
		++pl;
	size = s.length() - pl;
	vector<string> bit(size);
	for (auto& i : bit){
		switch (s[pl]){
		case '0': i = "0000"; break; case '4': i = "0100"; break; case '8': i = "1000"; break; case 'C': i = "1100"; break;
		case '1': i = "0001"; break; case '5': i = "0101"; break; case '9': i = "1001"; break; case 'D': i = "1101"; break;
		case '2': i = "0010"; break; case '6': i = "0110"; break; case 'A': i = "1010"; break; case 'E': i = "1110"; break;
		case '3': i = "0011"; break; case '7': i = "0111"; break; case 'B': i = "1011"; break; case 'F': i = "1111"; break;
		}
		++pl;
	}
	//Now, we have all bits of k
	BigNumber A(a);
	if (bit[size - 1][3] == '1')
		b = a;
	for (int i = 1; i < deg; i++){
		A *= A;
		A %= n;
		if (bit[size - 1 - i / 4][3 - i % 4] == '1'){
			b *= A;
			b %= n;
		}
	}
	return b;
}

BigNumber QuadraticSieve::LegendreSymbol(BigNumber& a, BigNumber& p){

	BigNumber r(modPow(a, (p - BigNumber::One()) / BigNumber::Two(), p));
	if (r > BigNumber::One())
		r -= p;
	return r;
}

std::pair<BigNumber, BigNumber> QuadraticSieve::doFactorization(){
	cout << "bit size of " << N << " is " << N.BitSize() << endl;
	double nLg = N.b_ln();
	double t1 = exp(sqrt(nLg) * sqrt(log(nLg)) * 0.5);
	Ipp32u t = ceil(t1);
	cout << "Factor base size chosen as "<<t << endl;

	Ipp32u len;
	if (N.BitSize() >= 200)
		len = t*ceill(sqrtl(sqrtl(t)));
	else if (N.BitSize() >= 50)
		len = t*ceill(sqrtl(t));  //TODO: optimize this bound
	else
		len = t*t;

	//Sieve of Eratosphenes
	std::vector<bool> arr(len + 1, true);
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

	Ipp32u counter = 1;
	Ipp32u l = len+1;
	Base = vector<BigNumber>(t+1);
	Base[0] = BigNumber::MinusOne();
	//This statemet are easy for parallel
	for (Ipp32u i = 2; i != l && counter <=t; ++i){
		if (arr[i]){
			BigNumber buf(i);
			if (LegendreSymbol(N, buf) == BigNumber::One()){      
				Base[counter] = buf;
				if (counter % 100000 == 0)
					cout << counter << endl;
				++counter;
			}
		}
	}

	if (counter - t > 1)
		cout << "We need more primes";

	return divisors;
}
