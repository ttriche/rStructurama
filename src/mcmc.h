#ifndef MCMC_H
#define MCMC_H

#include <ctime>
#include <fstream>
#include <string>
#include <vector>

using namespace std;



class MbRandom;
class McmcSamples;
class Model;
class NexusFile;
class Partition;
class Mcmc {

	public:
                            Mcmc(NexusFile &settings);
                            ~Mcmc(void);
                   double   getHeat(int index, double temperature) { return 1.0 / (1.0 + index * temperature); }
					 void   runChain(void);

	private:
					 void   attemptSwap(int **si);
	                  int   findColdChain(void);
				   string   formatTime(float n);
				     void   printToScreen(int n);
					 void   printSwapSummary(int **si);
					 void   sampleChain(int n);
				NexusFile   *settingsPtr;
	                  int   numCycles;
	                  int   printFrequency;
	                  int   sampleFrequency;
	                  int   numChains;
	               double   temperature;
				    Model   **models;
                  clock_t   startCpuTime;
				  clock_t   endCpuTime;
	          McmcSamples   *mcmcSamplesPtr;
			     MbRandom   *ranPtr;
	             ofstream   outputFile;
};

#endif
