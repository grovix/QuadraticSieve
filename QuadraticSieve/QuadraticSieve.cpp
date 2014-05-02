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
	double t1 = exp(0.5*sqrt(nLg) * sqrt(log(nLg))); //TODO: experiment with it
	Ipp32u t = ceil(t1);
	t /= 4;  //Optimize this
	//t = 4823;
	fbSize = t;

	cout << "First factor base size chosen as "<<fbSize << endl;

	Ipp32u len;
	if (N.BitSize() >= 200)
		len = t*ceill(sqrtl(sqrtl(t)));
	else if (N.BitSize() >= 50)
		len = t*ceill(sqrtl(t));  //TODO: optimize this bound
	else
		len = t*t;

	//Sieve of Eratosphenes
	std::vector<bool> arr(EratospheneSieve(len + 1));

	Ipp32u counter = 1;
	Ipp32u l = len-1;
	Base.push_back(BigNumber::MinusOne());
	Base.push_back(BigNumber::Two());
	//This statemet are easy for parallel
	for (Ipp32u i = 3; i != l && counter <t-1; ++i){
		if (arr[i]){
			BigNumber buf(i);
			if (LegendreSymbol(N, buf) == BigNumber::One()){      
				Base.push_back(buf);
				if (counter % 100000 == 0)
					cout << counter << endl;
				++counter;
			}
		}
	}

	cout << "Base size: " << Base.size() << endl;
	sieving();

	return divisors;
}

vector<pair<BigNumber, vector<Ipp32u>>> QuadraticSieve::sieving(){
	ofstream out("log.txt");
	float epsilon = 6;
	int prime_deg = 1; //Play with it

	Ipp32u decimal_size = (ceil((float)N.BitSize() / log2(10)));
	Ipp32u M;
	cout << "decimal length of number " << decimal_size << endl;
	//I'll optimize it later
	if (decimal_size <= 40)
		M = 100000; 
	else if (decimal_size <= 50)
		M = 200000;
	else if (decimal_size <= 55)
		M = 400000;
	else if (decimal_size <= 60)
		M = 500000;
	else if (decimal_size <= 67)
		M = 700000;
	else
		M = 1500000;

	//experiment
	BigNumber k = Base[Base.size() - 1];
	vector<Ipp32u> v;
	k.num2vec(v);
	M = v[0];
	cout << v.size() << endl;
	cout << " M " << M << endl;


	for (auto h = Base.begin(); h != Base.end(); ++h){
		if (*h > M){
			cout << h - Base.begin() << endl;
			break;
		}
	}

	//define pseudo random generator
	int ctxSize;
	int maxBitSize = ((N.BitSize()) / 2) - log2(M) +1;
	ippsPrimeGetSize(maxBitSize, &ctxSize);
	IppsPrimeState* pPrimeG = (IppsPrimeState*)(new Ipp8u[ctxSize]);
	ippsPrimeInit(maxBitSize, pPrimeG);

	ippsPRNGGetSize(&ctxSize);
	IppsPRNGState* pRand = (IppsPRNGState*)(new Ipp8u[ctxSize]);
	ippsPRNGInit(160, pRand);
	srand(time(NULL));
	ippsPRNGSetSeed(BN(BigNumber(rand())), pRand);

	//from this moment i can parallel my programm later

	//fill array of log(p[i])
	vector<float> prime_log(Base.size() - 1);

	auto b_end = Base.end();
	for (auto it = Base.begin() + 1; it != b_end; ++it){
		Ipp32u ind = it - Base.begin()-1;
		vector<Ipp32u> v;
		it->num2vec(v);
		prime_log[ind] = log(v[0]);
	}
	//result matrix
	vector<pair<BigNumber, vector<Ipp32u>>> result(Base.size() + 1);
	for (auto& i : result){
		i.second = std::vector<Ipp32u>(Base.size()+1);
	}

	Ipp32u counter = 0;
	Ipp32u bound = Base.size();
	BigNumber three(3);
	BigNumber four(4);
	BigNumber  D, h2, B2, primeMod, r1, r2, A2;

	BigNumber test;
	test = (N * BigNumber::Two()).b_sqrt() / M;
	cout << "approximate value of A " << test << endl;
	while (counter < bound){
		vector<float> sieve(2 * M + 1, 0);
		//generate coefficients
		ippsPrimeGen_BN(A, maxBitSize, nTrials, pPrimeG, ippsPRNGen, pRand);
		Ipp32u counter2 = 0;
		while (!(A.isPrime(nTrials) && LegendreSymbol(N,A) == BigNumber::One())){
			ippsPrimeGen_BN(A, maxBitSize, nTrials, pPrimeG, ippsPRNGen, pRand);
			counter2++;
			if (counter2 > 10000){ //optional
				maxBitSize += 1;
				int ctxSize;
				ippsPrimeGetSize(maxBitSize, &ctxSize);
				IppsPrimeState* pPrimeG = (IppsPrimeState*)(new Ipp8u[ctxSize]);
				ippsPrimeInit(maxBitSize, pPrimeG);

				ippsPRNGGetSize(&ctxSize);
				IppsPRNGState* pRand = (IppsPRNGState*)(new Ipp8u[ctxSize]);
				ippsPRNGInit(160, pRand);
				srand(time(NULL));
				ippsPRNGSetSeed(BN(BigNumber(rand())), pRand);

				counter2 = 0;
			}
		}
		cout << "Switch polynom" << endl;
		out << "Switch polynom" << endl;
		cout << "A " << A << endl;
		B = Tonelli_Shanks(N, A);
		cout << "B " << B << endl;
		C = (B*B - N) / A;
		//Q(x) = Ax^2 + 2*B*x +C
		//пока есть поиск корней только по простому модулю
		for (auto it = Base.begin() + 2; it != b_end-1; ++it){
			primeMod = BigNumber(*it);
			A2 = BigNumber::Two()*A;
			//cout << "primeModDeg " << primeModDeg << endl;
			BigNumber InvA;
			ippsModInv_BN(BN(A2 % primeMod), BN(primeMod), BN(InvA));
			D = Tonelli_Shanks((B*B - A*C), primeMod);
			r1 = InvA*(BigNumber::MinusOne()*BigNumber::Two()*B + BigNumber::Two()*D) % primeMod;
			r2 = InvA*(BigNumber::MinusOne()*BigNumber::Two()*B - BigNumber::Two()*D) % primeMod;
			Ipp32u ind = it - Base.begin() - 1;
			vector<Ipp32u> v, v1, v2;
			r1.num2vec(v);
			primeMod.num2vec(v1);
			r2.num2vec(v2);
			Ipp32s stepR1 = v[0];
			Ipp32u primePower = v1[0];
			Ipp32s stepR2 = v2[0];
			//cout << primePower << endl;
			if (primePower <= M){
				//sieving by first root
				while (stepR1 <= M){
					sieve[stepR1 + M] -= prime_log[ind];
					stepR1 += primePower;
				}
				stepR1 = v[0] + M - primePower;
				while (stepR1 >= 0){
					sieve[stepR1] -= prime_log[ind];
					stepR1 -= primePower;
				}
				//sieving by second root
				while (stepR2 <= M){
					sieve[stepR2 + M] -= prime_log[ind];
					stepR2 += primePower;
				}
				stepR2 = v2[0] + M - primePower;
				while (stepR2 >= 0){
					sieve[stepR2] -= prime_log[ind];
					stepR2 -= primePower;
				}
			}
		} //end for one prime from factor base
		auto s_end = sieve.end();

		for (auto&& it = sieve.begin(); it != sieve.end(); ++it){
			Ipp32s x = it - sieve.begin() - M;
			*it += Q(BigNumber(x)).b_abs().b_ln();
		}
		cout << "end initialize sieve" << endl;

		BigNumber R, Qx;
		for (auto it = sieve.begin(); it != s_end && counter <bound; ++it){
			    if (abs(*it) <= epsilon){
				//cout << "try" << endl;
				Ipp32s xg = it - sieve.begin() - M;
				BigNumber x(xg);
				Qx = Q(x);
				//cout << Qx << endl;

				BigNumber xt(Qx);
				if (Qx < BigNumber::Zero()){
					result[counter].second[0] = 1;
					Qx *= BigNumber::MinusOne();
				}
				//Qx %= N;
				//Trial division stage
				if (Qx > BigNumber::Zero()){
					for (auto jt = Base.begin() + 1; jt != Base.end(); ++jt){
						Ipp32u ind = jt - Base.begin();
						Ipp32u deg = 0;
						//Silverman point us to one improvement : R = x mod pi must bi equals one of two roots
						//I'll do it later
						while (Qx % *jt == BigNumber::Zero()){
							++deg;
							Qx /= *jt;
						}
						result[counter].second[ind] = deg;
					}
					if (Qx == BigNumber::One()){
						result[counter].first = (BigNumber::Two() * A * x + B)/ (BigNumber::Two() *D);
						cout << "counter ++ " <<counter<< endl;
						cout << "epsilon " << *it << endl;
						cout << xg << endl;
						out << "counter ++ " << counter << endl;
						out << "epsilon " << *it << endl;
						out << xg << endl;
						cout << result[counter].first << endl;
						++counter;
					}
					else{
						for (auto& kk : result[counter].second)  //I want to zero all degrees
							kk = 0;
					}
				}
			}
		}
	}

	return result;
}

BigNumber QuadraticSieve::Q(BigNumber& x){
	return A*x*x + 2 * B* x + C;
}

//Algorithm return x -> x^2 = a (mod p) or return false
BigNumber QuadraticSieve::Tonelli_Shanks(BigNumber& a, BigNumber& p){
	Ipp32u e = 0;
	BigNumber q = p - BigNumber::One();
	while (q.IsEven()){
		++e;
		q /= 2;
	}
	//1.Find generator
	BigNumber n;
	int size;
	int numSize = 10;  //What size "n" should be?

	ippsPRNGGetSize(&size);
	IppsPRNGState* pPrng = (IppsPRNGState*)(new Ipp8u[size]);
	ippsPRNGInit(160, pPrng);

	ippsPRNGSetSeed(BN(BigNumber(rand())), pPrng);

	ippsPRNGen_BN(BN(n), numSize, pPrng);
	while (LegendreSymbol(n, p) != BigNumber::MinusOne())
		ippsPRNGen_BN(BN(n), numSize, pPrng);

	BigNumber z(modPow(n, q, p));
	//2.Initialize
	BigNumber two(2);
	BigNumber t;
	BigNumber y(z);
	BigNumber r(e);
	BigNumber x(modPow(a, (q - BigNumber::One()) / BigNumber::Two(), p));
	BigNumber b(a*(x*x %p) % p);
	x = a*x %p;
	while (true){
		if (b % p == BigNumber::One())
			return x;
		//3.Find exponent
		BigNumber m(1);
		while (modPow(b, two.b_power(m), p) != BigNumber::One())
			m += 1;
		if (m == r)
			return BigNumber::Zero();
		//4.Reduce exponent
		t = modPow(y, two.b_power(r - m - BigNumber::One()), p);
		y = t*t;
		y %= p;
		r = m;
		x *= t;
		x %= p;
		b *= y;
		b %= p;
	}

}
