#include "QuadraticSieve.h"
#include "wiedemann.h"
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
	BigNumber F(a);
	if (bit[size - 1][3] == '1')
		b = a;
	for (int i = 1; i < deg; i++){
		F *= F;
		F %= n;
		if (bit[size - 1 - i / 4][3 - i % 4] == '1'){
			b *= F;
			b %= n;
		}
	}
	return b;
}

BigNumber QuadraticSieve::LegendreSymbol(const BigNumber& a, const BigNumber& p){

	BigNumber r(modPow(a, (p - BigNumber::One()) / BigNumber::Two(), p));
	if (r > BigNumber::One())
		r -= p;
	return r;
}

std::pair<BigNumber, BigNumber> QuadraticSieve::doFactorization(){

	std::pair<BigNumber, BigNumber> divisors;
	cout << "bit size of " << N << " is " << N.BitSize() << endl;
	double nLg = N.b_ln();
	double t1 = exp(0.5*sqrt(nLg) * sqrt(log(nLg))); //TODO: experiment with it
	Ipp32u t = ceil(t1);
	t /= 8;  //Optimize this
	if (N.BitSize() > 150)
		t /= 4;
	if (N.BitSize() > 200)
		t /= 4;

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
	fbSize = Base.size();
	vector<pair<BigNumber, vector<Ipp32u>>> smooth(sieving());

	cout << "Sieving stage completed" << endl;
	SparseMatrix SM(getSparseMatrix(smooth), fbSize + 1);
	Wiedemann calc(SM);

	vector<bool> w(calc.getSolution());
	while (SM.isZero(w))
		w = calc.getSolution();

	bool isCompleted = false;
	while (!isCompleted){
		BigNumber x(1);
		for (uInt i = 0; i < fbSize + 1; ++i){
			if (w[i]){
				x *= smooth[i].first;
				//x %= N;
			}
		}
		x %= N;
		cout <<"x "<< x << endl;
		BigNumber y(1);
		vector<Ipp32u> L(fbSize,0);
		for (Ipp32u j = 0; j < fbSize; ++j){
			for (Ipp32u i = 0; i < fbSize + 1; ++i){
				if (w[i])
					L[j] += smooth[i].second[j];
			}
			L[j] /= 2;	
		}
		for (Ipp32u j = 0; j < fbSize; ++j){
			y *= Base[j].b_power(L[j]);
		}
		y %= N;
		BigNumber d(N);
		IppStatus f1 = ippsGcd_BN(BN(x-y), BN(N), BN(d));
		cout << "f1 " << f1 << endl;
		if (d != N && d != BigNumber::One() && d != BigNumber::Zero()){
			cout << "Completed !!! " << endl;
			isCompleted = true;
			divisors.first = d;
			divisors.second = N / d;
		}
		else{
			cout << "d " << d << endl;
			w = calc.getSolution();
			while (SM.isZero(w))
				w = calc.getSolution();
		}
	}
	return divisors;
}

vector<pair<BigNumber, vector<Ipp32u>>> QuadraticSieve::sieving(){
	float epsilon = 6;

	Ipp32u decimal_size = (ceil((float)N.BitSize() / log2(10)));
	Ipp32s M;
	cout << "decimal length of number " << decimal_size << endl;
	//I'll optimize it

	BigNumber k = Base[fbSize - 1];
	vector<Ipp32u> v;
	k.num2vec(v);
	M = v[0];
	M *= 3;
	M /= 2;
	cout << " M " << M << endl;

	//define pseudo random generator
	int ctxSize;
	int maxBitSize = ((N.BitSize()) / 2) - log2(M)+1;
	maxBitSize /= 2;
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
	vector<pair<BigNumber, vector<Ipp32u>>> result(fbSize + 1);
	for (auto& i : result){
		i.second = std::vector<Ipp32u>(fbSize+1);
	}

	Ipp32u counter = 0;
	Ipp32u bound = fbSize+1;
	BigNumber three(3);
	BigNumber four(4);
	BigNumber  D, h2, B2, primeMod, r1, r2, A2;

	BigNumber test;
	BigNumber q;
	float* sieve = new float[2 * M + 1];
	for (int i = 0; i < 2 * M + 1; ++i)
		sieve[i] = 0;
	while (counter < bound){
		//generate coefficients
		ippsPrimeGen_BN(q, maxBitSize, nTrials, pPrimeG, ippsPRNGen, pRand);
		while (!(q.isPrime(nTrials) && LegendreSymbol(N,q) == BigNumber::One())){
			ippsPrimeGen_BN(q, maxBitSize, nTrials, pPrimeG, ippsPRNGen, pRand);
		}
		cout << "counter = " << counter << endl;
		//cout << "Switch polynom" << endl;
		A = q*q;

		//now use Hensel�s Lemma
		BigNumber B1 = Tonelli_Shanks(N, q);
		BigNumber invF(q);
		ippsModInv_BN(BN(BigNumber::Two()*B1 %q), BN(q), BN(invF));

		BigNumber t = ((N - B1*B1) / q)*invF % q;
		B = (B1 + q*t) % A;
		C = (B*B - N) / A;
		BigNumber Q1(N);
		ippsModInv_BN(BN(q), BN(N), BN(Q1));
		//Q(x) = Ax^2 + 2*B*x +C
		//���� ���� ����� ������ ������ �� �������� ������
		for (auto it = Base.begin() + 2; it != b_end ; ++it){
			primeMod = BigNumber(*it);
			A2 = BigNumber::Two()*A;
			BigNumber InvA(primeMod);
			ippsModInv_BN(BN(A2 % primeMod), BN(primeMod), BN(InvA));
			D = Tonelli_Shanks(N, primeMod);
			r1 = InvA*(BigNumber::MinusOne()*BigNumber::Two()*B + BigNumber::Two()*D) % primeMod;
			r2 = InvA*(BigNumber::MinusOne()*BigNumber::Two()*B - BigNumber::Two()*D) % primeMod;

			Ipp32u ind = it - Base.begin() - 1;
			vector<Ipp32u> v, v1, v2;
			r1.num2vec(v);
			primeMod.num2vec(v1);
			r2.num2vec(v2);
			Ipp32s stepR1 = v[0];
			Ipp32s primePower = v1[0];
			Ipp32s stepR2 = v2[0];
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
		//cout << "all roots calculated"<<endl;
		//clock_t start = clock();
		for (Ipp32s i = 0; i <2*M+1; ++i){
			Ipp32s x = i - M;
			sieve[i] += Q(BigNumber(x)).b_abs().b_ln();
		}
		
		BigNumber R, Qx;
		for (Ipp32s i = 0; i <2 * M + 1 & counter < bound; ++i){
			if (sieve[i] <= epsilon){
				Ipp32s xg = i - M;
				BigNumber x(xg);
				Qx = Q(x);

				BigNumber xt(Qx);
				if (Qx < BigNumber::Zero()){
					result[counter].second[0] = 1;
					Qx *= BigNumber::MinusOne();
				}
				//Trial division stage
				if (Qx > BigNumber::Zero()){
					for (auto jt = Base.begin() + 1; jt != Base.end(); ++jt){
						Ipp32u ind = jt - Base.begin();
						Ipp32u deg = 0;
						//Silverman point us to one improvement : R = x mod pi must be equals one of two roots
						//I'll do it later
						while (Qx % *jt == BigNumber::Zero()){
							++deg;
							Qx /= *jt;
						}
						result[counter].second[ind] = deg;
					}
					if (Qx == BigNumber::One()){
						result[counter].first = (A*x + B) * Q1; //???
						++counter;
					}
					else{
						for (auto& kk : result[counter].second)  //I want to zero all degrees
							kk = 0;
					}
				}
			}
		}
		for (int i = 0; i < 2 * M + 1; ++i)
			sieve[i] = 0;
	}
	return std::move(result);
}

BigNumber QuadraticSieve::Q(const BigNumber& x){
	return std::move(A*x*x + 2 * B* x + C);
}

//Algorithm return x -> x^2 = a (mod p) or return false
BigNumber QuadraticSieve::Tonelli_Shanks(const BigNumber& a,const BigNumber& p){
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
			m += BigNumber::One();
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

vector<unsigned int> QuadraticSieve::getSparseMatrix(vPair& v){
	unsigned int len = v.size();
	vector<unsigned int> r;
	for (unsigned int i = 0; i < len; ++i){
		r.push_back(0);
		const auto ind = r.end() - r.begin() - 1;
		for (unsigned int j = 0; j < len; ++j){
			if (v[j].second[i] % 2 != 0){
				r.push_back(j);
				r[ind]++;
			}
		}
	}
	return std::move(r);
}