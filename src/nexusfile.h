#ifndef NEXUSFILE_H
#define NEXUSFILE_H

#include <istream>
#include <vector>
#include <map>

using namespace std;



class CompInterface;
class IoManager;
class Locus;
class McmcSamples;
class NexusBlock;
class Observations;
class UserInterface;
class NexusFile {

	public:
                            NexusFile(void);
                            ~NexusFile(void);
							
                     bool   deleteSample(void);
            CompInterface   *getCompInterfacePtr(void) { return ci; }
                IoManager   *getChainOutPtr(void) { return chainOut; }
			  McmcSamples   *getMcmcSamplesPtr(void) { return mcmcSamples; }
				IoManager   *getMeanPartOutPtr(void) { return meanPartOut; }
				IoManager   *getNumPopsOutPtr(void) { return numPopsOut; }
				IoManager   *getTogethernessOutPtr(void) { return togethernessOut; }
                     void   interpretCmd(istream &c, UserInterface *ui);
                     bool   makeNewSample(void);
             Observations   *samplePtr(void) { return sample; }
                     void   showRecodedObservations(void);
                     void   setAdmixtureCompInterfacePtr(CompInterface *p) { ci = p; }
                     void   setCompInterfacePtr(CompInterface *p) { ci = p; }
					 
					 void   setAdmixture(bool tf) { admixture = tf; }
				     void   setAdmixtureDist(string s) { admixtureDist = s; }
				     void   setAdmixtureParm1(double x) { admixtureParm1 = x; }
				     void   setAdmixtureParm2(double x) { admixtureParm2 = x; }
				     void   setAdmixtureConcentrationDist(string s) { admConcentrationPriorDist = s; }
				     void   setAdmixtureConcentrationParm1(double x) { admConcentrationParm1 = x; }
				     void   setAdmixtureConcentrationParm2(double x) { admConcentrationParm2 = x; }
					 void   setBurnInMarginalLike(int x) { burnInMarginalLike = x; }
                     void   setBurnInMeanPart(int x) { burnInMeanPart = x; }
                     void   setBurnInNumPops(int x) { burnInNumPops = x; }
                     void   setBurnInTogetherness(int x) { burnInTogehterness = x; }
				     void   setConcentrationDist(string s) { concentrationDist = s; }
				     void   setConcentrationParm1(double x) { concentrationParm1 = x; }
				     void   setConcentrationParm2(double x) { concentrationParm2 = x; }
                     void   setNgen(int x) { ngen = x; }
                     void   setNumChains(int x) { numChains = x; }
                     void   setNumIndividuals(int x) { numIndividuals = x; }
                     void   setNumLoci(int x) { numLoci = x; }
                     void   setNumPops(int x) { numPops = x; }
                     void   setPrintfreq(int x) { printfreq = x; }
                     void   setSamplefreq(int x) { samplefreq = x; }
                     void   setSampleRead(bool tf) { sampleRead = tf; }
					 void   setSaveToFileMeanPart(bool x) { saveToFileMeanPart = x; }
					 void   setSaveToFileNumPops(bool x) { saveToFileNumPops = x; }
					 void   setSaveToFileTogetherness(bool x) { saveToFileTogetherness = x; }
                     void   setSumMeanpart(bool tf) { calcMeanPart = tf; }
                     void   setTemperature(double x) { temperature = x; }
			         void   setWarnUser(bool tf) { warnUser = tf; }
					 
					 bool   getAdmixture(void) { return admixture; }
				   string   getAdmixtureDist(void) { return admixtureDist; }
				   double   getAdmixtureParm1(void) { return admixtureParm1; }
				   double   getAdmixtureParm2(void) { return admixtureParm2; }
				   string   getAdmixtureConcentrationDist(void) { return admConcentrationPriorDist; }
				   double   getAdmixtureConcentrationParm1(void) { return admConcentrationParm1; }
				   double   getAdmixtureConcentrationParm2(void) { return admConcentrationParm2; }
					  int   getBurnInMarginalLike(void) { return burnInMarginalLike; }
                      int   getBurnInMeanPart(void) { return burnInMeanPart; }
                      int   getBurnInTogetherness(void) { return burnInTogehterness; }
                      int   getBurnInNumPops(void) { return burnInNumPops; }
				   string   getConcentrationDist(void) { return concentrationDist; }
				   double   getConcentrationParm1(void) { return concentrationParm1; }
				   double   getConcentrationParm2(void) { return concentrationParm2; }
                      int   getNgen(void) { return ngen; }
                      int   getNumChains(void) { return numChains; }
                      int   getNumIndividuals(void) { return numIndividuals; }
                      int   getNumLoci(void) { return numLoci; }
                      int   getNumPops(void) { return numPops; }
                      int   getPrintfreq(void) { return printfreq; }
                   double   getTemperature(void) { return temperature; }
                      int   getSamplefreq(void) { return samplefreq; }
                     bool   getSampleRead(void) { return sampleRead; }
					 bool   getSaveToFileMeanPart(void) { return saveToFileMeanPart; }
					 bool   getSaveToFileNumPops(void) { return saveToFileNumPops; }
					 bool   getSaveToFileTogetherness(void) { return saveToFileTogetherness; }
                     bool   getSumMeanpart(void) { return calcMeanPart; }
					 bool   getWarnUser(void) { return warnUser; }

	private:
				   string   admConcentrationPriorDist;
				   double   admConcentrationParm1;
				   double   admConcentrationParm2;
				     bool   admixture;
				   string   admixtureDist;
				   double   admixtureParm1;
				   double   admixtureParm2;
					  int   burnInMarginalLike;
                      int   burnInMeanPart;
					  int   burnInNumPops;
					  int   burnInTogehterness;
                     bool   calcMeanPart;
                IoManager   *chainOut;
            CompInterface   *ci;
				   string   concentrationDist;
				   double   concentrationParm1;
				   double   concentrationParm2;
			  McmcSamples   *mcmcSamples;
				IoManager   *meanPartOut;
	 vector<NexusBlock *>   nexusBlocks;
                      int   ngen;
                      int   numChains;
                      int   numIndividuals;
					  int   numLoci;
                      int   numPops;
				IoManager   *numPopsOut;
                      int   printfreq;
                     void   readNexusID(istream &c);
             Observations   *sample;
                      int   samplefreq;
                     bool   sampleRead;
					 bool   saveToFileMeanPart;
					 bool   saveToFileNumPops;
					 bool   saveToFileTogetherness;
                   double   temperature;
				IoManager   *togethernessOut;
				     bool   warnUser;
};



class NexusCommand : public string {

	public:
                            NexusCommand (istream &);
};



class NexusBlock {

	public:
                     bool   readSetting(istream &c, string &key, string &val);
                     void   readKeyValPairs(istream &c, map<string, string> &keyValPairs);
};  



class NexusDataBlock : public NexusBlock {

	public:
                            NexusDataBlock(NexusFile &parent, istream &c, UserInterface *ui); 

	private:
                     bool   dimensionsFound;
                     void   doDimensions(NexusFile &parent, istream &c, UserInterface *ui);
                     bool   doInfo(NexusFile &parent, istream &c, UserInterface *ui);
};



class NexusFile;
class UserInterface;
class NexusStructuramaBlock : public NexusBlock {

	public:
                            NexusStructuramaBlock(NexusFile &parent, istream &c, UserInterface *ui);
                     void   doAcknowledge(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doCalcMarginal(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doCitation(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doExecute(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doSet(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doShowdata(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doHelp(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doMcmc(NexusFile &parent, istream &str, UserInterface *ui); 
                     void   doModel(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doReadSamples(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doShowMeanPart(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doShowNumPops(NexusFile &parent, istream &str, UserInterface *ui);
                     void   doShowTogetherness(NexusFile &parent, istream &str, UserInterface *ui);

	private:
					  int   interpretDistribution(string &qs, string &dist, double &val1, double &val2);
             static const   string parmNames[];
             static const   vector<string> parms; 
               static int   numParms();
};





#endif
