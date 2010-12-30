#ifndef POPPARTITION_H
#define POPPARTITION_H

#include <vector>

using namespace std;


class Franchise;
class MbBitfield;
class Restaurant;
class PopPartition {

	public:
                            PopPartition(Restaurant *rp);
                            PopPartition(Franchise *fp);
							PopPartition(PopPartition &p);
							PopPartition(vector<PopPartition *> &partList, double &score, vector<double> &var, bool &cancelMeanPartitionProcess);
							PopPartition(vector<int> &rgf);
                            ~PopPartition(void);
		     PopPartition   &operator=(PopPartition &p);
				     bool   operator==(const PopPartition &P);
					 void   deleteElement(int elemToDelete);
					  int   distance(PopPartition *p);
                      int   getDegree(void) { return degree; }
					  int   getElement(int i) { return part[i]; }
					  int   getNumElements(void) { return numElements; }
					  int   *getPart(void) { return part; }
			  MbBitfield*   getSubset(int i) { return subsets[i]; }
					 void   print(void);
					 void   setElement(int i, int j) { part[i] = j; }
					 void   setRgf(void);
					 void   setDegree(void);
					 void   setSubsets(void);
				   double   ssDistance(vector<PopPartition *> &partList);

	private:
	                 void   allocatePartition(void);
					 void   deleteSubsets(void);
					  int   flip(int x);
					 void   retreiveDiv(int whichDiv, int *div);
				   double   ssDistance(vector<PopPartition *> &partList, int elemToDelete);
	                  int   numElements;
					  int   degree;
	                  int   *part;
	 vector<MbBitfield *>   subsets;
};

#endif
