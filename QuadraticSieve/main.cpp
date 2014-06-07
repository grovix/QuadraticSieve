//////////////////////////////////////////////////
/* Integer factorization using quadratic sieve
*  with Intel Compiler and Intel IPP
*
*  Grigoriy Trofimyuk, 2014
*/////////////////////////////////////////////////

#include <fstream>
#include <ipp.h>
#include <stdio.h>
#include "Factorization.h"
#include <time.h>
#include "wiedemann.h"
#include <array>
#include <thread>
vector<BigNumber> BigNumber::decPowers;

void testFill(int**a, int id){
	for (int i = 0; i < 4; i++)
		a[id][i] = id;
}

void threadTest(){

}
int main(){
	fstream in("input.txt");
	std::string num;
	in >> num;

	clock_t start = clock();

	Factorization comp = Factorization(num.c_str());
	BigNumber test(num.c_str());

	//инициализируем вектор степеней десятки
	Ipp32u deg = ceil((test.BitSize() - 1) / log2f(10));
	++deg;
	BigNumber ten(10);
	BigNumber tmp(1);
	for (unsigned i = 0; i <= deg; ++i){
		BigNumber::decPowers.push_back(tmp);
		tmp *= ten;
	}

	std::map<BigNumber, Ipp32u> factor = comp.getFactor();
	cout << "Factorization completed! " << endl<<"time = "<<clock() - start << endl;
	for (auto& i : factor){
		cout << i.first;
		if (i.second > 1)
			cout << "^" << i.second;
		cout << endl;
	}

	//QuadraticSieve q(test);
	//q.doFactorization();

	cout << "Complete " << clock() - start << endl;
	ofstream out("TS_log.txt");
	out << "Complete " << clock() - start << endl;
	system("pause");
	return 0;
}