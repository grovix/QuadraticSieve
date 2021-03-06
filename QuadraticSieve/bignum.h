#if !defined _BIGNUMBER_H_
#define _BIGNUMBER_H_
#include <ipp.h>
#include <ippcp.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <time.h>


using namespace std;
class BigNumber
{
public:
	BigNumber(Ipp32u value = 0);
	BigNumber(Ipp32s value);
	BigNumber(const IppsBigNumState* pBN);
	BigNumber(const Ipp32u* pData, int length = 1, IppsBigNumSGN sgn = IppsBigNumPOS);
	BigNumber(const BigNumber& bn);
	BigNumber(const char *s);
	virtual ~BigNumber();
	// set value
	void Set(const Ipp32u* pData, int length = 1, IppsBigNumSGN sgn = IppsBigNumPOS);
	// conversion to IppsBigNumState
	friend IppsBigNumState* BN(const BigNumber& bn) { return bn.m_pBN; }
	operator IppsBigNumState* () const { return m_pBN; }
	// some useful constatns
	static const BigNumber& Zero();
	static const BigNumber& One();
	static const BigNumber& Two();
	static const BigNumber& MinusOne();
	// arithmetic operators probably need
	BigNumber& operator = (const BigNumber& bn);
	BigNumber& operator += (const BigNumber& bn);
	BigNumber& operator -= (const BigNumber& bn);
	BigNumber& operator *= (Ipp32u n);
	BigNumber& operator *= (const BigNumber& bn);
	BigNumber& operator /= (const BigNumber& bn);
	BigNumber& operator %= (const BigNumber& bn);
	friend BigNumber operator + (const BigNumber& a, const BigNumber& b);
	friend BigNumber operator - (const BigNumber& a, const BigNumber& b);
	friend BigNumber operator * (const BigNumber& a, const BigNumber& b);
	friend BigNumber operator * (const BigNumber& a, Ipp32u);
	friend BigNumber operator % (const BigNumber& a, const BigNumber& b);
	friend BigNumber operator / (const BigNumber& a, const BigNumber& b);
	// modulo arithmetic
	BigNumber Modulo(const BigNumber& a) const;
	BigNumber ModAdd(const BigNumber& a, const BigNumber& b) const;
	BigNumber ModSub(const BigNumber& a, const BigNumber& b) const;
	BigNumber ModMul(const BigNumber& a, const BigNumber& b) const;
	BigNumber InverseAdd(const BigNumber& a) const;
	BigNumber InverseMul(const BigNumber& a) const;
	// comparisons
	friend bool operator < (const BigNumber& a, const BigNumber& b);
	friend bool operator > (const BigNumber& a, const BigNumber& b);
	friend bool operator == (const BigNumber& a, const BigNumber& b);
	friend bool operator != (const BigNumber& a, const BigNumber& b);
	friend bool operator <= (const BigNumber& a, const BigNumber& b) { return !(a>b); }
	friend bool operator >= (const BigNumber& a, const BigNumber& b) { return !(a<b); }
	// easy tests
	bool IsOdd() const;
	bool IsEven() const { return !IsOdd(); }
	// size of BigNumber
	int MSB() const;
	int LSB() const;
	int BitSize() const { return MSB() + 1; }
	int DwordSize() const { return (BitSize() + 31) >> 5; }
	// conversion and output
	void num2hex(string& s) const; // convert to hex string
	void num2vec(vector<Ipp32u>& v) const; // convert to 32-bit word vector
	friend ostream& operator << (ostream& os, const BigNumber& a);
	//my modifiņations
	bool isPrime(int nTrials);
	float b_ln();
	float b_log2();
	BigNumber b_sqrt();
	BigNumber b_gcd(const BigNumber& a);
	BigNumber b_power(const BigNumber& e);
	BigNumber b_abs();
	static vector<BigNumber> decPowers;
	//
protected:
	bool create(const Ipp32u* pData, int length, IppsBigNumSGN sgn = IppsBigNumPOS);
	int compare(const BigNumber&) const;
	IppsBigNumState* m_pBN;
};

inline int Bit(const vector<Ipp32u>& v, int n)
{
	return 0 != (v[n >> 5] & (1 << (n & 0x1F)));
}
// convert bit size into 32-bit words
#define BITSIZE_WORD(n) ((((n)+31)>>5))
#endif // _BIGNUMBER_H_

vector<bool> EratospheneSieve(Ipp32u bound);