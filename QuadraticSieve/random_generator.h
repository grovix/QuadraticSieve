#if !defined _RANDOM_GENERATOR_H
#define _RANDOM_GENERATOR_H

#include <ippcp.h>
#include <iostream>
class RandGen{
public:
	RandGen(){
		status = ippsPRNGGetSize(&PRNGStateSize);
		if (status != ippStsNoErr)
			throw("Error in allocation context of random generator");
		status = ippsPRNGInit(seedBits, PRNGStateContext);
		if (status != ippStsNoErr)
			throw("Error in initialization PRNG");

	}
	int PRNGStateSize;
	IppsPRNGState* PRNGStateContext = (IppsPRNGState*)(new Ipp8u[PRNGStateSize]);
	IppStatus status;
	int seedBits = 17; //magic number =)
private:

};

#endif