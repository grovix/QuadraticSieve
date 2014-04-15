#include <fstream>
#include "bignum.h"
#include <ipp.h>
#include <stdio.h>
#include "Factorization.h"

int main(){
	int n;
	std::cout << "Choose input" << std::endl <<
		"1-read from file" << std::endl <<
		"2-read from keyboard" << std::endl;
	cin >> n;
	fstream in("input.txt");
	string num;
	switch (n){
	case 1:
		in >> num;
	}
	system("pause");
	return 0;
}