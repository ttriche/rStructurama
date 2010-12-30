#ifndef STIRLING_H
#define STIRLING_H

#include "integer.h"


class Stirling {

	public:
	                        Stirling(int n);
							~Stirling(void);
				   double   getLnStirling(int i) { return lnStirling[i]; }
					 void   calculateProbK(double alpha);
					 void   calculateLnProbK(double alpha);
				   double   getProbK(int i) { return priorOnK[i]; }
				   double   getLnProbK(int i) { return priorOnK[i]; }

	private:
	                  int   num;
				   double   *lnStirling;
				   double   *priorOnK;
};

#endif



