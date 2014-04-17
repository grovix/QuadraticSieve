#include <fstream>
#include "bignum.h"
#include <ipp.h>
#include <stdio.h>
#include "Factorization.h"

int main(){

	fstream in("input.txt");
	std::string num;
	in >> num;
	Factorization comp = Factorization(num.c_str());

	comp.checkPrime();
	system("pause");
	return 0;
}