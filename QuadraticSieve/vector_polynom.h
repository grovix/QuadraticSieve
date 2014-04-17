#if !defined _VECTOR_POLYNOM_H
#define _VECTOR_POLYNOM_H

#include "polynom.h"

class VectorPolynom :public Polynom{
public:
	VectorPolynom();
	VectorPolynom(int degree);
	VectorPolynom(vector<bool> poly);
	friend VectorPolynom operator + (const VectorPolynom& a, const VectorPolynom& b);
	bool subs(int degree);
private:
	int p_degree;
	vector<bool> polynom;
};



#endif