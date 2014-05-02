#include "bignum.h"
//////////////////////////////////////////////////////////////////////
//
// BigNumber
//
//////////////////////////////////////////////////////////////////////
BigNumber::~BigNumber()
{
	delete[](Ipp8u*)m_pBN;
}
bool BigNumber::create(const Ipp32u* pData, int length, IppsBigNumSGN sgn)
{
	int size;
	ippsBigNumGetSize(length, &size);
	m_pBN = (IppsBigNumState*)(new Ipp8u[size]);
	if (!m_pBN)
		return false;
	ippsBigNumInit(length, m_pBN);
	if (pData)
		ippsSet_BN(sgn, length, pData, m_pBN);
	return true;
}
// constructors
//
BigNumber::BigNumber(Ipp32u value)
{
	create(&value, 1, IppsBigNumPOS);
}
BigNumber::BigNumber(Ipp32s value)
{
	Ipp32s avalue = abs(value);
	create((Ipp32u*)&avalue, 1, (value<0) ? IppsBigNumNEG :
		IppsBigNumPOS);
}
BigNumber::BigNumber(const IppsBigNumState* pBN)
{
	int bnLen;
	ippsGetSize_BN(pBN, &bnLen);
	Ipp32u* bnData = new Ipp32u[bnLen];
	IppsBigNumSGN sgn;
	ippsGet_BN(&sgn, &bnLen, bnData, pBN);
	//
	create(bnData, bnLen, sgn);
	//
	delete bnData;
}
BigNumber::BigNumber(const Ipp32u* pData, int length, IppsBigNumSGN sgn)
{
	create(pData, length, sgn);
}
static
char HexDigitList[] = "0123456789ABCDEF";
BigNumber::BigNumber(const char* s)
{
	bool neg = '-' == s[0];
	if (neg) s++;
	bool hex = ('0' == s[0]) && (('x' == s[1]) || ('X' == s[1]));
	int dataLen;
	Ipp32u base;
	if (hex) {
		s += 2;
		base = 0x10;
		dataLen = (strlen(s) + 7) / 8;
	}
	else {
		base = 10;
		dataLen = (strlen(s) + 9) / 10;
	}
	create(0, dataLen);
	*(this) = Zero();
	while (*s) {
		char tmp[2] = { s[0], 0 };
		Ipp32u digit = strcspn(HexDigitList, tmp);
		*this = (*this) * base + BigNumber(digit);
		s++;
	}
	if (neg)
		(*this) = Zero() - (*this);
}
BigNumber::BigNumber(const BigNumber& bn)
{
	IppsBigNumSGN sgn;
	int length;
	ippsGetSize_BN(bn.m_pBN, &length);
	Ipp32u* pData = new Ipp32u[length];
	ippsGet_BN(&sgn, &length, pData, bn.m_pBN);
	//
	create(pData, length, sgn);
	//
	delete[] pData;
}
// set value
//
void BigNumber::Set(const Ipp32u* pData, int length, IppsBigNumSGN sgn)
{
	ippsSet_BN(sgn, length, pData, BN(*this));
}
// constants
//
const BigNumber& BigNumber::Zero()
{
	static const BigNumber zero(0);
	return zero;
}
const BigNumber& BigNumber::One()
{
	static const BigNumber one(1);
	return one;
}
const BigNumber& BigNumber::Two()
{
	static const BigNumber two(2);
	return two;
}

const BigNumber& BigNumber::MinusOne()
{
	static const BigNumber minus_one(-1);
	return minus_one;
}
// arithmetic operators
//
BigNumber& BigNumber::operator =(const BigNumber& bn)
{
	if (this != &bn) { // prevent self copy
		int length;
		ippsGetSize_BN(bn.m_pBN, &length);
		Ipp32u* pData = new Ipp32u[length];
		IppsBigNumSGN sgn;
		ippsGet_BN(&sgn, &length, pData, bn.m_pBN);
		//
		delete (Ipp8u*)m_pBN;
		create(pData, length, sgn);
		//
		delete pData;
	}
	return *this;
}
BigNumber& BigNumber::operator += (const BigNumber& bn)
{
	int aLength;
	ippsGetSize_BN(BN(*this), &aLength);
	int bLength;
	ippsGetSize_BN(BN(bn), &bLength);
	int rLength = IPP_MAX(aLength, bLength) + 1;
	BigNumber result(0, rLength);
	ippsAdd_BN(BN(*this), BN(bn), BN(result));
	*this = result;
	return *this;
}
BigNumber& BigNumber::operator -= (const BigNumber& bn)
{
	int aLength;
	ippsGetSize_BN(BN(*this), &aLength);
	int bLength;
	ippsGetSize_BN(BN(bn), &bLength);
	int rLength = IPP_MAX(aLength, bLength);
	BigNumber result(0, rLength);
	ippsSub_BN(BN(*this), BN(bn), BN(result));
	*this = result;
	return *this;
}
BigNumber& BigNumber::operator *= (const BigNumber& bn)
{
	int aLength;
	ippsGetSize_BN(BN(*this), &aLength);
	int bLength;
	ippsGetSize_BN(BN(bn), &bLength);
	int rLength = aLength + bLength;
	BigNumber result(0, rLength);
	ippsMul_BN(BN(*this), BN(bn), BN(result));
	*this = result;
	return *this;
}
BigNumber& BigNumber::operator *= (Ipp32u n)
{
	int aLength;
	ippsGetSize_BN(BN(*this), &aLength);
	BigNumber bn(n);
	BigNumber result(0, aLength + 1);
	ippsMul_BN(BN(*this), BN(bn), BN(result));
	*this = result;
	return *this;
}
BigNumber& BigNumber::operator %= (const BigNumber& bn)
{
	BigNumber remainder(bn);
	ippsMod_BN(BN(*this), BN(bn), BN(remainder));
	*this = remainder;
	return *this;
}
BigNumber& BigNumber::operator /= (const BigNumber& bn)
{
	BigNumber quotient(*this);
	BigNumber remainder(bn);
	ippsDiv_BN(BN(*this), BN(bn), BN(quotient), BN(remainder));
	*this = quotient;
	return *this;
}
BigNumber operator + (const BigNumber& a, const BigNumber &b)
{
	BigNumber r(a);
	return r += b;
}
BigNumber operator - (const BigNumber& a, const BigNumber &b)
{
	BigNumber r(a);
	return r -= b;
}
BigNumber operator * (const BigNumber& a, const BigNumber &b)
{
	BigNumber r(a);
	return r *= b;
}
BigNumber operator * (const BigNumber& a, Ipp32u n)
{
	BigNumber r(a);
	return r *= n;
}
BigNumber operator / (const BigNumber& a, const BigNumber &b)
{
	BigNumber q(a);
	return q /= b;
}
BigNumber operator % (const BigNumber& a, const BigNumber &b)
{
	BigNumber r(b);
	ippsMod_BN(BN(a), BN(b), BN(r));
	return r;
}
// modulo arithmetic
//
BigNumber BigNumber::Modulo(const BigNumber& a) const
{
	return a % *this;
}
BigNumber BigNumber::InverseAdd(const BigNumber& a) const
{
	BigNumber t = Modulo(a);
	if (t == BigNumber::Zero())
		return t;
	else
		return *this - t;
}
BigNumber BigNumber::InverseMul(const BigNumber& a) const
{
	BigNumber r(*this);
	ippsModInv_BN(BN(a), BN(*this), BN(r));
	return r;
}
BigNumber BigNumber::ModAdd(const BigNumber& a, const BigNumber &b) const
{
	BigNumber r = this->Modulo(a + b);
	return r;
}
BigNumber BigNumber::ModSub(const BigNumber& a, const BigNumber &b) const
{
	BigNumber r = this->Modulo(a + this->InverseAdd(b));
	return r;
}
BigNumber BigNumber::ModMul(const BigNumber& a, const BigNumber &b) const
{
	BigNumber r = this->Modulo(a*b);
	return r;
}
// comparison
//
int BigNumber::compare(const BigNumber &bn) const
{
	Ipp32u result;
	BigNumber tmp = *this - bn;
	ippsCmpZero_BN(BN(tmp), &result);
	return (result == IS_ZERO) ? 0 : (result == GREATER_THAN_ZERO) ? 1 : -1;
}
bool operator < (const BigNumber &a, const BigNumber &b)
{
	return a.compare(b) < 0;
}
bool operator > (const BigNumber &a, const BigNumber &b)
{
	return a.compare(b) > 0;
}
bool operator == (const BigNumber &a, const BigNumber &b)
{
	return 0 == a.compare(b);
}
bool operator != (const BigNumber &a, const BigNumber &b)
{
	return 0 != a.compare(b);
}
// easy tests
//
bool BigNumber::IsOdd() const
{
	vector<Ipp32u> v;
	num2vec(v);
	return v[0] & 1;
}
// size of BigNumber
//
int BigNumber::LSB() const
{
	if (*this == BigNumber::Zero())
		return 0;
	vector<Ipp32u> v;
	num2vec(v);
	int lsb = 0;
	vector<Ipp32u>::iterator i;
	for (i = v.begin(); i != v.end(); i++) {
		Ipp32u x = *i;
		if (0 == x)
			lsb += 32;
		else {
			while (0 == (x & 1)) {
				lsb++;
				x >>= 1;
			}
			break;
		}
	}
	return lsb;
}
int BigNumber::MSB() const
{
	if (*this == BigNumber::Zero())
		return 0;
	vector<Ipp32u> v;
	num2vec(v);
	int msb = v.size() * 32 - 1;
	vector<Ipp32u>::reverse_iterator i;
	for (i = v.rbegin(); i != v.rend(); i++) {
		Ipp32u x = *i;
		if (0 == x)
			msb -= 32;
		else {
			while (!(x & 0x80000000)) {
				msb--;
				x <<= 1;
			}
			break;
		}
	}
	return msb;
}
int Bit(const vector<Ipp32u>& v, int n)
{
	return 0 != (v[n >> 5] & (1 << (n & 0x1F)));
}
// conversions and output
//
void BigNumber::num2vec(vector<Ipp32u>& v) const
{
	int length;
	ippsGetSize_BN(BN(*this), &length);
	Ipp32u* pData = new Ipp32u[length];
	IppsBigNumSGN sgn;
	ippsGet_BN(&sgn, &length, pData, BN(*this));
	//
	for (int n = 0; n<length; n++)
		v.push_back(pData[n]);
	//
	delete pData;
}
void BigNumber::num2hex(string& s) const
{
	int length;
	ippsGetSize_BN(BN(*this), &length);
	Ipp32u* pData = new Ipp32u[length];
	IppsBigNumSGN sgn;
	ippsGet_BN(&sgn, &length, pData, BN(*this));
	s.append(1, (sgn == IppsBigNumNEG) ? '-' : ' ');
	s.append(1, '0');
	s.append(1, 'x');
	for (int n = length; n>0; n--) {
		Ipp32u x = pData[n - 1];
		for (int nd = 8; nd>0; nd--) {
			char c = HexDigitList[(x >> (nd - 1) * 4) & 0xF];
			s.append(1, c);
		}
	}
	delete pData;
}
ostream& operator << (ostream &os, const BigNumber& a) {
	string s;
	a.num2hex(s);
	os << s.c_str();
	return os;
}

bool BigNumber::isPrime(int nTrials){

	//define Number-bit Prime Generator
	bool res;
	int ctxSize;
	int bitSize = this->BitSize();
	ippsPrimeGetSize(bitSize, &ctxSize);
	IppsPrimeState* pPrimeG = (IppsPrimeState*)(new Ipp8u[ctxSize]);
	ippsPrimeInit(bitSize, pPrimeG);

	// define Pseudo Random Generator (default settings)
	ippsPRNGGetSize(&ctxSize);
	IppsPRNGState* pRand = (IppsPRNGState*)(new Ipp8u[ctxSize]);
	ippsPRNGInit(160, pRand);

	Ipp32u result;
	ippsPrimeTest_BN(BN(*this), nTrials, &result, pPrimeG, ippsPRNGen, pRand);
	if (result == IPP_IS_PRIME)
		res = true;
	else
		res = false;

	delete[](Ipp8u*)pRand;
	delete[](Ipp8u*)pPrimeG;

	return res;
}

float BigNumber::b_ln(){
	vector<Ipp32u> v;
	if (*this < 20000){
		this->num2vec(v);
		return log((float)v[0]);
	}

	int n = this->BitSize() - 1;
	float deg = ceilf(n / log2f(10)); //evaluation of decimal points
	deg -= 4;
	n = (int)deg;
	char* dTen = new char[n + 1];
	dTen[0] = '1';
	for (int i = 1; i < n; ++i){
		dTen[i] = '0';
	}
	dTen[n] = '\0';

	(*this / dTen).num2vec(v);
	delete [] dTen;

	return (log((float)v[0]) + ((float)n-1)*log(10)); //log(1000) = 3*log(10)
}

BigNumber BigNumber::b_sqrt(){ //Newton's method

	vector<Ipp32u> v;
	if (*this < 10000){
		this->num2vec(v);
		return std::move(BigNumber((Ipp32u) sqrt((float)v[0])));
	}

	int n = this->BitSize() - 1;
	float deg = ceilf(n / log2f(10)); //evaluation of decimal points
	deg -= 4;
	n = (int)deg;
	n /= 2;
	char* dTen = new char[n + 1];
	dTen[0] = '1';
	for (int i = 1; i < n; ++i){
		dTen[i] = '0';
	}
	dTen[n] = '\0';

	BigNumber x(dTen);
	BigNumber xn(0);

	delete[] dTen;

	while (true){
		xn = (x + *this / x);
		xn /= Two();
		if (xn < x){
			if (x - xn < Two())
				break;
		}
		else if ((xn - x) < Two())
			break;
		
		x = xn;
	}
	return std::move(xn);
}

BigNumber BigNumber::b_gcd(const BigNumber& a){
	BigNumber res;
	ippsGcd_BN(BN(*this), BN(a), BN(res));
	return res;
}

vector<bool> EratospheneSieve(Ipp32u bound){
	vector<bool> arr(bound,true);
	auto& a_begin = arr.begin();
	auto& a_end = arr.end();
	Ipp32u ind, p;
	arr[0] = false;
	arr[1] = false;
	auto& ibound = a_begin + sqrt(a_end - a_begin);
	for (auto& i = arr.begin() + 2; i <= ibound; ++i){
		if (*i == true){
			p = i - a_begin;
			ind = p * p;
			while (ind < (a_end - a_begin)){
				arr[ind] = false;
				ind += p;
			}
		}
	}
	return std::move(arr);
}

BigNumber BigNumber::b_power(const BigNumber& e){
	BigNumber p(e);
	BigNumber res(1);
	while(p > BigNumber::Zero()){
		res *= *this;
		p -= 1;
	}
	return res;
}

BigNumber BigNumber::b_abs(){
	if (*this < BigNumber::Zero())
		return BigNumber::MinusOne()* (*this);
	else
		return *this;
}