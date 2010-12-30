#ifndef SAMPLES_H
#define SAMPLES_H

#include <vector>
#include <string>

using namespace std;

struct PatronId {

				      int   individualId;
				      int   locusId;
				      int   alleleId;
					 bool   isIndividual;
				   double   latitude;
				   double   longitude;
};

class NexusFile;
class PopPartition;
class Sample {

	public:
                            Sample(void);
                            ~Sample(void);
					 void   addAdmixtureVariancePrior(string s) { admixtureVariancePrior = s; }
					 void   addAdmixtureVariance(float x) { admixtureVariance = x; }
					 void   addCycle(int x) { cycle = x; }
					 void   addLnLike(float x) { lnLike = x; }
					 void   addDppAlphaPopulation(float x) { dppAlphaPopulation = x; }
					 void   addDppAlphaPopulationPrior(string s) { dppAlphaPopulationPrior = s; }
					 void   addDppAlphaAdmixture(float x) { dppAlphaAdmixture = x; }
					 void   addDppAlphaAdmixturePrior(string s) { dppAlphaAdmixturePrior = s; }
					 void   addPartition(PopPartition *p) { partition = p; }
				    float   getLnLike(void) { return lnLike; }
			 PopPartition   *getPartition(void) { return partition; }
			       string   getSampleString(void);

	private:
	                float   admixtureVariance;
				   string   admixtureVariancePrior;
	                  int   cycle;
					float   dppAlphaPopulation;
				   string   dppAlphaPopulationPrior;
			        float   dppAlphaAdmixture;
				   string   dppAlphaAdmixturePrior;
				    float   lnLike;
			 PopPartition   *partition;
};

class Model;
class McmcSamples {

	public:
                            McmcSamples(void);
                            ~McmcSamples(void);
					 void   addHeader(string s) { header.push_back( s ); }
					 void   addPatronId(PatronId *p) { patronIds.push_back( p ); }
					 void   addSample(Sample *s) { samples.push_back( s ); }
					 void   calcMarginalLike(NexusFile &settings);
					 void   calcMeanPart(NexusFile &settings);
					 void   calcTogetherness(NexusFile &settings);
					 void   calcNumPops(NexusFile &settings);
					 void   deleteSamples(void);
					 bool   getGpsInfoPresent(void) { return gpsInfoPresent; }
				   string   getHeaderString(Model *mp);
				    float   *getJointProbs(void) { return jointProbs; }
				   double   getLatitude(int i);
				   double   getLongitude(int i);
					  int   getNumGpsCoordinates(void);
					  int   getNumIndividuals(void);
					  int   getNumSamples(void) { return samples.size(); }
					  int   getNumPatrons(void);
					  int   getNumMixturesForPatronOnMeanPartition(int i);
					  int   getDegreeOfMeanPartition(void);
					  int   getSmallestPatronId(void);
					  int   getLargestPatronId(void);
					float   *getMeanPartitionMixtureProbs(void);
					 void   getMixtureProbsForMeanPartition(int indIdx, float *probs, int numMixtures);
				 PatronId   *getPatronId(int i) { return patronIds[i]; }
			 PopPartition   *getMeanPartition(void) { return meanPartition; }
			         bool   hasGps(PatronId *p);
					 bool   isEmpty(void);
					 bool   readFile(IoManager &mngr, NexusFile &settings);
					 void   setPatronIds(Model *mp);

	private:
	               double   bestMeanPartScore;
				     bool   gpsInfoPresent;
	       vector<string>   header;
	                 void   interpretInd(string h, int &a, int &b, int &c, bool &d, double &e, double &f);
			        float   *jointProbs;
			 PopPartition   *meanPartition;
			        float   *meanPartitionMixtureProbs;
				      int   numSamplesMeanPartIsBasedOn;
	   vector<PatronId *>   patronIds;
	     vector<Sample *>   samples;
		   vector<double>   var;
};

#endif
