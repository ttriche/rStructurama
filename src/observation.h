#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <string>
#include <list>
#include <vector>

using namespace std;


class MbBitfield;
class Observations;
class Individual {

	public:
                            Individual(Observations *op, int nl);
                            ~Individual(void);
                     void   allocInfo(int x);
                     void   freeInfo(void);
					  int   getInfoAt(int l, int a);
				   double   getLatitude(void) { return latitude; }
				   double   getLongitude(void) { return longitude; }
                   string   getName(void) { return name; }
                      int   getNumAllelesAtLocus(int l) { return numAlleles[l]; }
                      int   getNumLoci(void) { return numLoci; }
                      int   getNumUniqueAllelesAtLocus(int i);
					  int   getNumUniqueAlleles(void);
                      int   getRawInfoAt(int l, int a) { return rawInfo[l][a]; }
					 bool   hasGpsCoordinates(void);
					 void   printRecoded(void);
					 void   setInfo(int l, int a, int s);
					 void   setLatitude(double x) { latitude = x; }
					 void   setLongitude(double x) { longitude = x; }
                     void   setName(string s) { name = s; }
                     void   setNumLoci(int x) { numLoci = x; }
                     void   setRawInfo(int loc, int a[2], int na) { rawInfo[loc][0]=a[0]; rawInfo[loc][1]=a[1]; numAlleles[loc]=na; }

	private:
	         Observations   *obsPtr;
                   string   name;
                      int   numLoci;
                      int   *numAlleles;
                      int   **rawInfo;
                     bool   isInfoAllocated;
				   double   latitude;
				   double   longitude;
			   MbBitfield   *info[2];
};



class Observations {

	public:
                            Observations(int ni, int nl);
                            ~Observations(void);
                     bool   getIsDiploid(int i) { return isDiploid[i]; }
                   string   getName(int i) { return individuals[i]->getName(); }
                      int   getNumIndividuals(void) { return numIndividuals; }
                      int   getNumUniqueAllelesAtLocus(int i) { return numAllelesAtLocus[i]; }
					  int   getNumUniqueAlleles(void) { return numUniqueAlleles; }
                      int   getNumLoci(void) { return numLoci; }
					  int   getPositionForLocus(int l) { return positionForLocus[l]; }
               Individual   *getPtrToIndividual(int i) { return individuals[i]; }
			          int   getTotalNumAlleleCopies(void) { return totalNumAlleleCopies; }
                     void   print(void);
                     void   printRecoded(void);
                     void   recodeInfo(void);
                     void   setNumIndividuals(int x) { numIndividuals = x; }
                     void   setNumLoci(int x) { numLoci = x; }
	
	private:
	                  int   numIndividuals;
     vector<Individual *>   individuals;
                list<int>   getListOfAllelesAtLocus(int l);
	                  int   numLoci;
					  int   numUniqueAlleles;
	                  int   *numAllelesAtLocus;
	                  int   *positionForLocus;
	                 bool   *isDiploid;
					  int   totalNumAlleleCopies;

};



#endif
