#if !defined _POLYNOM_H
#define _POLYNOM_H

#include <ippcp.h>
#include <iostream>
#include <vector>
#include <iterator>
using namespace std;

class Polynom
{
public:
	friend Polynom operator + (const Polynom& a, const Polynom& b);
	bool subs(int degree);
};


#endif