#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>
#include "MbVector.h"


using namespace std;


class Franchise;
class Individual;
class MbRandom;
class Model;
class NexusFile;
class Observations;
class Restaurant;
class Table;

class Patron {

	public:
                            Patron(Observations *op, int pIdx, int iIdx, int nde);
							~Patron(void);
					  int   getIndex(void) { return patronIdx; }
					  int   getIndividual(void) { return individualIdx; }
					  int   getDataOffset(void) { return offSet; }
					 int*   getDataPtr(void) { return data; }
				   double   getLatitude(void) { return gps[0]; }
				   double   getLongitude(void) { return gps[1]; }
				   double   getLnLikelihoodalone(void) { return lnLikelihoodAlone; }
					  int   getNumDataElements(void) { return numDataElements; }
					 bool   getPatronIsIndividual(void) { return patronIsIndividual; }
				   Table*   getTable(void) { return tablePtr; }
				     bool   hasGpsCoordinates(void);
				     void   setGps(double x, double y) { gps[0] = x; gps[1] = y; }
					 void   setIndex(int x) { patronIdx = x; }
					 void   setIndividual(int x) { individualIdx = x; }
			         void   setLnLikelihoodAlone(double x) { lnLikelihoodAlone = x; }
					 void   setTable(Table *t) { tablePtr = t; }

			 virtual void   print(void)=0;
			 virtual void   printData(void)=0;
			  virtual int   getLocus(void)=0;
			  virtual int   getAllele(void)=0;
		   virtual string   getPatronStr(void)=0;
			 virtual void   setData(Individual *ip)=0;
			 
	protected:
			 Observations   *observationsPtr;
	                Table   *tablePtr;
					  int   patronIdx;
					  int   individualIdx;
					  int   numDataElements;
					  int   *data;
					  int   offSet;
					 bool   patronIsIndividual;
				   double   lnLikelihoodAlone;
				   double   gps[2];
};

class PatronAllele : public Patron {

	public:
                            PatronAllele(Observations *op, int pIdx, int iIdx, int lIdx, int aIdx, int nde);
					 void   print(void);
					 void   printData(void);
					 void   setData(Individual *ip);
					  int   getLocus(void) { return locusIdx; }
					  int   getAllele(void) { return alleleIdx; }
				   string   getPatronStr(void);
	private:
					  int   locusIdx;
					  int   alleleIdx;
};

class PatronIndividual : public Patron {

	public:
                            PatronIndividual(Observations *op, int pIdx, int nl, int nde);
					 void   print(void);
					 void   printData(void);
					 void   setData(Individual *ip);
					  int   getLocus(void) { return 0; }
					  int   getAllele(void) { return 0; }
				   string   getPatronStr(void);

	private:
					  int   numLoci;
};

class Menu;
class Table {

	public:
                            Table(Restaurant *rsp, MbRandom *rp, Observations *op, int nde, double *ft);
							~Table(void);
					 void   addObservationsForPatron(Patron *p);
					 void   subtractObservationsForPatron(Patron *p);
					 void   addPatronToTable(Patron *p);
					  int   getNumPatrons(void) { return seatedPatrons.size(); }
				  Patron*   getPatron(int i) { return seatedPatrons[i]; }
					 int*   getTableAlleleCountsPtr(void) { return &alleleCounts[0]; }
				    Menu*   getMenuItem(void) { return menuPtr; }
				   double   lnLikelihood(Patron *p);
				   double   lnLikelihood(void);
					 void   print(void);
					 void   printData(void);
					 void   removePatronFromTable(Patron *p);
					 void   setMenuItem(Menu *m) { menuPtr = m; }

	private:
	           Restaurant   *restaurantPtr;
				 MbRandom   *ranPtr;
			 Observations   *observationsPtr;
		 vector<Patron *>   seatedPatrons;
		             Menu   *menuPtr;
					  int   numDataElements;
		              int   *alleleCounts;
				   double   *lnFactorial;
};

class Franchise;
class Menu {

	public:
                            Menu(Franchise *fp, MbRandom *rp, Observations *op, int nde, double *ft);
							~Menu(void);
				     void   addObservationsForPatron(Patron *p);
					 void   addTableToMenu(Table *t);
					 void   deleteTableFromMenu(Table *t);
					  int   getNumTables(void) { return seatedTables.size(); }
				   Table*   getTable(int i) { return seatedTables[i]; }
				   double   lnLikelihood(void);
				   double   lnLikelihood(Patron *p);
					 void   print(void);
					 void   printData(void);
					 void   subtractObservationsForPatron(Patron *p);

	private:
	            Franchise   *franchisePtr;
				 MbRandom   *ranPtr;
			 Observations   *observationsPtr;
		  vector<Table *>   seatedTables;
					  int   numDataElements;
		              int   *alleleCounts;
				   double   *lnFactorial;
};

class Restaurant {

	public:
                            Restaurant(Model *mp, MbRandom *rp, Observations *op, NexusFile *sp, bool dpp, bool adm, int wi, double cp);
							~Restaurant(void);
					 void   addTableToRestaurant(Table *tp);
	                 void   deleteTable(Table *tp);
				   double   getAdmixtureVarianceParm(void) { return admixtureVarianceParm; }
				   double   getConcentrationParm(void) { return concentrationParm; }
					  int   getNumTables(void) { return tables.size(); }
				   Table*   getTable(int i) { return tables[i]; }
				      int   getNumPatronsInRestaurant(void) { return numPatrons; }
				  Patron*   getRestaurantPatron(int i) { return patrons[i]; }
				  double*   getFactorialPtr(void) { return lnFactorial; }
				   double   lnLikelihood(void);
				   double   lnPriorProbability(void);
					 void   normalizeVector(vector<double> &v);
					 void   print(void);
					 void   printData(void);
					  int   multinomialRv(vector<double> &v);
					 void   setConcentrationParm(double x) { concentrationParm = x; }
					 void   updateAdmixtureProportions(void);
					 void   updateAdmixtureVarianceParm(void);
					 void   updateConcParm(void);
					 void   updateConcParm(double a, double b);
					 void   updateSeatingForDpp(void);
					 void   updateSeatingForFixedTables(void);
					 
	private:
				      int   getIndexForTable(Table *tp);
	                 void   initializeAdmixtureProportions(void);
	                 void   initializeFactorials(void);
					 void   initializePatrons(void);
					 void   initializeTables(void);
				   double   lnLikelihoodNewTable(Patron *p);
				   double   concentrationParm;
				     bool   assumingAdmixture;
					 bool   assumingDpp;
					  int   individualsInRestaurant;
	                Model   *modelPtr; 
				 MbRandom   *ranPtr;
			 Observations   *observationsPtr;
					  int   numPatrons;
		           Patron   **patrons;
			    NexusFile   *settingsPtr;
		  vector<Table *>   tables;
				   double   *lnFactorial;
			MbVector<int>   *popCounts;
		 MbVector<double>   *popProbs;
		           double   admixtureVarianceParm;
				   double   admixtureLowerLimit;
				   double   admixtureUpperLimit;
};

class Stirling;
class Franchise {

	public:
                            Franchise(Model *mp, MbRandom *rp, Observations *op, NexusFile *sp);
							~Franchise(void);
				      int   getNumPatronsInFranchise(void);
					  int   numRestaurantTables(void);
					  int   getNumRestaurants(void) { return restaurants.size(); }
			  Restaurant*   getRestaurant(int i) { return restaurants[i]; }
			          int   getIndexForMenuItem(Menu *mp);
				   double   getMenuAlpha(void) { return menuAlpha; }
					  int   getNumMenuItems(void) { return menus.size(); }
					  int   getNumTablesInFranchise(void);
				   double   getRestaurantAlpha(void) { return restaurantAlpha; }
				   double   lnLikelihood(void);
				   double   lnPriorProbability(void);
					 void   print(void);
					 void   printMenuData(void);
				   double   safeExponentiation(double lnR);
					 void   updateConcParm(void);
					 void   updateSeating(void);

	private:
	                 void   deleteMenuItem(Menu *mToDelete);
				   double   lnLikelihoodNewTable(Patron *p, vector<double> &pv);
				   double   sumLogProbs(vector<double> &v);
	                Model   *modelPtr; 
				 MbRandom   *ranPtr;
			 Observations   *observationsPtr;
			    NexusFile   *settingsPtr;
	 vector<Restaurant *>   restaurants;
		   vector<Menu *>   menus;
				   double   restaurantAlpha;
				   double   menuAlpha;
				 Stirling   *stirlingNums;
};

class PopPartition;
class Model {

	public:
                            Model(int idx, NexusFile &settings);
							~Model(void);
					  int   getNumElements(void) { return numElements; }
					  int   getNumElementsPerRestaurant(void) { return numElementsPerRestaurant; }
					  int   getNumPopulations(void) { return numPopulations; }
					  int   getCurrentNumberOfPopulations(void);
				   double   getAlpha1(void) { return alpha1; }
				   double   getAlpha2(void) { return alpha2; }
			   Franchise*   getFranchisePtr(void) { return franchise; }
				      int   getModelId(void) { return modelId; }
					  int   getModelIndex(void) { return modelIndex; }
			PopPartition*   getPartition(void);
			  Restaurant*   getRestaurantPtr(void) { return restaurant; }
			       string   getConcentrationParm1Model(void) { return concentrationParm1Model; }
				   string   getConcentrationParm2Model(void) { return concentrationParm2Model; }
				   string   getAdmixtureModel(void) { return admixtureModel; }
			       double   lnLikelihood(void);
				   double   lnPriorProbability(void);
					 void   setAlpha1(double x) { alpha1 = x; }
					 void   setModelIndex(int x) { modelIndex = x; }
				   double   getAdmixtureVarianceParm(void);
				     bool   isNumPopsRandom(void);
					 void   print(void);
					 void   updateSeating(void);
					 void   updateAdmixtureProportions(void);
					 void   updateAdmixtureVarianceParm(void);
					 void   updateConcParm(void);

	private:
	               double   calcAlphaFromExpectedNumberOfTables(int n, double expT);
				   double   expNumTables(int n, double a);
				   double   expNumTables(int n, double alpha, double beta);
				 MbRandom   *ranPtr;
			 Observations   *observationsPtr;
	                  int   modelIndex;
					  int   modelId;
					  int   numElements;
					  int   numElementsPerRestaurant;
			          int   numPopulations;
				   double   alpha1;
				   double   conc1Shape;
				   double   conc1Scale;
				   double   conc1Lambda;
				   double   alpha2;
				   double   conc2Shape;
				   double   conc2Scale;
				   double   conc2Lambda;
				   double   alpha3;
				   double   admShape;
				   double   admScale;
				   double   admLambda;
				   double   admUniform[2];
	               string   concentrationParm1Model;
	               string   concentrationParm2Model;
	               string   admixtureModel;
		       Restaurant   *restaurant;
			    Franchise   *franchise;
};




#endif
