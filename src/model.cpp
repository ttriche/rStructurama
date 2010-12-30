#include "MbRandom.h"
#include "MbVector.h"
#include "model.h"
#include "nexusfile.h"
#include "observation.h"
#include "partition.h"
#include "st.h"
#include "stirling.h"
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <numeric>
#include <algorithm>

#undef NO_DATA

using namespace std;




#pragma mark Patron (Base Class)

Patron::Patron(Observations *op, int pIdx, int iIdx, int nde) {

	observationsPtr = op;
	tablePtr        = NULL;
	patronIdx       = pIdx;
	individualIdx   = iIdx;
	numDataElements = nde;
	gps[0]          = 1000.0;
	gps[1]          = 1000.0;
	if (op->getPtrToIndividual(iIdx)->hasGpsCoordinates() == true)
		setGps( op->getPtrToIndividual(iIdx)->getLatitude(), op->getPtrToIndividual(iIdx)->getLongitude() );
	data = new int[numDataElements];
	for (int i=0; i<numDataElements; i++)
		data[i] = 0;
}

Patron::~Patron(void) {

	delete [] data;
}

bool Patron::hasGpsCoordinates(void) {

	if (gps[0] > -90.0 && gps[0] < 90.0 && gps[1] > -180.0 && gps[1] < 180.0)
		return true;
	return false;
}

#pragma mark Patron (Allele)

PatronAllele::PatronAllele(Observations *op, int pIdx, int iIdx, int lIdx, int aIdx, int nde) : Patron(op, pIdx, iIdx, nde) {

	locusIdx           = lIdx;
	alleleIdx          = aIdx;
	patronIsIndividual = false;
}

void PatronAllele::print(void) {

	cout << "(" << getIndividual() << "," << locusIdx << "," << alleleIdx << ")";
}

string PatronAllele::getPatronStr(void) {

	char tempStr[100];
	sprintf(tempStr, "(%d,%d,%d)", getIndividual(), locusIdx, alleleIdx);
	string str(tempStr);
	return str;
}

void PatronAllele::printData(void) {
	
	for (int i=0; i<locusIdx; i++)
		{
		for (int j=0; j<observationsPtr->getNumUniqueAllelesAtLocus(i); j++)
			cout << "    ";
		}
	int *dataPtr = &data[0];
	for (int i=0; i<numDataElements; i++)
		{
		cout << (*dataPtr);
		dataPtr++;
		}
	cout << endl;
}

void PatronAllele::setData(Individual *ip) {

	data[ ip->getInfoAt(locusIdx, alleleIdx) ]++;
	offSet = 0;
	for (int l=0; l<locusIdx; l++)
		offSet += observationsPtr->getNumUniqueAllelesAtLocus(l);
}

#pragma mark Patron (Individual)

PatronIndividual::PatronIndividual(Observations *op, int pIdx, int nl, int nde) : Patron(op, pIdx, pIdx, nde) {

	numLoci            = nl;
	patronIsIndividual = true;
}

string PatronIndividual::getPatronStr(void) {

	char tempStr[100];
	sprintf(tempStr, "(%d)", getIndividual());
	string str(tempStr);
	return str;
}

void PatronIndividual::print(void) {

	cout << getIndividual();
}

void PatronIndividual::printData(void) {

	int *dataPtr = &data[0];
	for (int i=0; i<numDataElements; i++)
		{
		cout << (*dataPtr);
		dataPtr++;
		}
	cout << endl;
}

void PatronIndividual::setData(Individual *ip) {

	int *dataPtr = &data[0];
	for (int l=0; l<numLoci; l++)
		{
		for (int a=0; a<ip->getNumAllelesAtLocus(l); a++)
			{
			int alleleId = ip->getInfoAt(l, a);
			if (alleleId >= 0)
				dataPtr[ alleleId ]++;
			}
		dataPtr += ip->getNumUniqueAllelesAtLocus(l);
		}
	offSet = 0;
}

#pragma mark Table

Table::Table(Restaurant *rsp, MbRandom *rp, Observations *op, int nde, double *ft) {

	restaurantPtr   = rsp;
	ranPtr          = rp;
	observationsPtr = op;
	numDataElements = nde;
	lnFactorial     = ft;
	menuPtr         = NULL;
	alleleCounts = new int[numDataElements];
	for (int i=0; i<numDataElements; i++)
		alleleCounts[i] = 0;
}

Table::~Table(void) {

	delete [] alleleCounts;
	seatedPatrons.clear();
}

void Table::addPatronToTable(Patron *p) {

	// add the patron's pointer to the vector of patron pointers for this table
	seatedPatrons.push_back( p );
	p->setTable(this);
	
	// add the patron's observations to this table's pooled observations
	addObservationsForPatron( p );
	
	// add the observations to the franchise menu, if available
	if (menuPtr != NULL)
		menuPtr->addObservationsForPatron( p );
}

void Table::addObservationsForPatron(Patron *p) {

	int *a = &alleleCounts[ p->getDataOffset() ];
	int *b = p->getDataPtr();
	for (int i=0; i<p->getNumDataElements(); i++)
		a[i] += b[i];
}

double Table::lnLikelihood(Patron *p) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnP = 0.0;
	int *y = &alleleCounts[ p->getDataOffset() ];
	int *x = p->getDataPtr();
	
	if ( p->getPatronIsIndividual() == true )
		{
		// calculate likelihood for all loci: no admixture
		for (int l=0; l<observationsPtr->getNumLoci(); l++)
			{
			int sum1 = 0, sum2 = 0;
			int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
			for (int j=0; j<k; j++)
				{
				lnP += ( lnFactorial[(*x) + (*y)] - lnFactorial[(*y)] );
				sum1 += (*y);
				sum2 += (*x) + (*y);
				x++;
				y++;
				}
			lnP += ( lnFactorial[sum1+k] - lnFactorial[sum2+k] );
			}
		}
	else
		{
		// calculate likelihood for one locus only: admixture
		int k = p->getNumDataElements();
		int y_kl = 0, y_klj = 0;
		for (int j=0; j<k; j++)
			{
			if (x[j] == 1)
				y_klj = y[j];
			y_kl += y[j];
			}
		lnP = log( y_klj + 1.0 ) - log( y_kl + k );
		}
	return lnP;
}

double Table::lnLikelihood(void) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnL = 0.0;
	int *x = &alleleCounts[0];
	for (int l=0; l<observationsPtr->getNumLoci(); l++)
		{
		int sum = 0;
		int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
		for (int j=0; j<k; j++)
			{
			lnL += lnFactorial[(*x)];
			sum += (*x);
			x++;
			}
		lnL += ( lnFactorial[k] - lnFactorial[sum] );
		}	
	return lnL;
}

void Table::removePatronFromTable(Patron *p) {

	// remove the patron's pointer from the vector of patron pointers for this table
	for (vector<Patron *>::iterator e=seatedPatrons.begin(); e != seatedPatrons.end(); e++)
		{
		if ( (*e) == p )
			{
			seatedPatrons.erase( e );
			break;
			}
		}
	p->setTable(NULL);

	// subtract the patron's observations from the table's
	subtractObservationsForPatron( p );
	
	// subtract the observations from the franchise table, if available
	if (menuPtr != NULL)
		menuPtr->subtractObservationsForPatron( p );
}

void Table::subtractObservationsForPatron(Patron *p) {

	int *a = &alleleCounts[ p->getDataOffset() ];
	int *b = p->getDataPtr();
	for (int i=0; i<p->getNumDataElements(); i++)
		a[i] -= b[i];
}

void Table::print(void) {

	for (vector<Patron *>::iterator p=seatedPatrons.begin(); p != seatedPatrons.end(); p++)
		{
		(*p)->print();
		cout << " ";
		}
	cout << endl;
}

void Table::printData(void) {

	for (vector<Patron *>::iterator p=seatedPatrons.begin(); p != seatedPatrons.end(); p++)
		{
		cout << setw(8) << (*p)->getIndividual() << ": ";
		(*p)->printData();
		}
	cout << "     sum: ";
	int *x = &alleleCounts[0];
	for (int l=0; l<observationsPtr->getNumLoci(); l++)
		{
		int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
		cout << "(";
		for (int j=0; j<k; j++)
			{
			cout << x[j];
			if (j+1 != k)
				cout << ",";
			}
		x += k;
		cout << ") ";
		}	
	cout << endl;
}

#pragma mark Menu

Menu::Menu(Franchise *fp, MbRandom *rp, Observations *op, int nde, double *ft) {

	franchisePtr    = fp;
	ranPtr          = rp;
	observationsPtr = op;
	numDataElements = nde;
	lnFactorial     = ft;
	alleleCounts = new int[numDataElements];
	for (int i=0; i<numDataElements; i++)
		alleleCounts[i] = 0;
}

Menu::~Menu(void) {

	delete [] alleleCounts;
}

void Menu::addTableToMenu(Table *t) {

	// add the restaurant table to the franchise table
	seatedTables.push_back( t );
	
	// add the restaurant allele counts to the allele counts for this franchise table
	int *a = &alleleCounts[0];
	int *b = t->getTableAlleleCountsPtr();
	for (int i=0; i<numDataElements; i++)
		a[i] += b[i];
}

void Menu::addObservationsForPatron(Patron *p) {

	int *a = &alleleCounts[ p->getDataOffset() ];
	int *b = p->getDataPtr();
	for (int i=0; i<p->getNumDataElements(); i++)
		a[i] += b[i];
}

void Menu::deleteTableFromMenu(Table *t) {

	// remove the restaurant table from this franchise menu
	for (vector<Table *>::iterator e=seatedTables.begin(); e != seatedTables.end(); e++)
		{
		if ( (*e) == t )
			{
			seatedTables.erase( e );
			break;
			}
		}
		
	// subtract the restaurant allele counts from the allele counts for this franchise menu
	int *a = &alleleCounts[0];
	int *b = t->getTableAlleleCountsPtr();
	for (int i=0; i<numDataElements; i++)
		a[i] -= b[i];
}

double Menu::lnLikelihood(void) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnL = 0.0;
	int *x = &alleleCounts[0];
	for (int l=0; l<observationsPtr->getNumLoci(); l++)
		{
		int sum = 0;
		int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
		for (int j=0; j<k; j++)
			{
			lnL += lnFactorial[(*x)];
			sum += (*x);
			x++;
			}
		lnL += ( lnFactorial[k] - lnFactorial[sum] );
		}	
	return lnL;
}

double Menu::lnLikelihood(Patron *p) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnP = 0.0;
	int *y = &alleleCounts[ p->getDataOffset() ];
	int *x = p->getDataPtr();
	
	if ( p->getPatronIsIndividual() == true )
		{
		// calculate likelihood for all loci
		for (int l=0; l<observationsPtr->getNumLoci(); l++)
			{
			int sum1 = 0, sum2 = 0;
			int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
			for (int j=0; j<k; j++)
				{
				lnP += ( lnFactorial[(*x) + (*y)] - lnFactorial[(*y)] );
				sum1 += (*y);
				sum2 += (*x) + (*y);
				x++;
				y++;
				}
			lnP += ( lnFactorial[sum1+k] - lnFactorial[sum2+k] );
			}
		}
	else
		{
		// calculate likelihood for one locus only
		int sum1 = 0, sum2 = 0;
		int k = p->getNumDataElements();
		for (int j=0; j<k; j++)
			{
			lnP += ( lnFactorial[(*x) + (*y)] - lnFactorial[(*y)] );
			sum1 += (*y);
			sum2 += (*x) + (*y);
			x++;
			y++;
			}
		lnP += ( lnFactorial[sum1+k] - lnFactorial[sum2+k] );
		}
	return lnP;
}

void Menu::print(void) {

	stprint("Tables seated at menu item %d: ", this );
	for (vector<Table *>::iterator t=seatedTables.begin(); t != seatedTables.end(); t++)
		cout << (*t) << " ";
	cout << endl;
}

void Menu::printData(void) {

	stprint("Data for menu item %d:", this );
	int sum = 0;
	for (int i=0; i<numDataElements; i++)
		{
		sum += alleleCounts[i];
		cout << alleleCounts[i] << " ";
		}
	cout << "(" << sum << ")";
}

void Menu::subtractObservationsForPatron(Patron *p) {

	int *a = &alleleCounts[ p->getDataOffset() ];
	int *b = p->getDataPtr();
	for (int i=0; i<p->getNumDataElements(); i++)
		a[i] -= b[i];
}

#pragma mark Restaurant

Restaurant::Restaurant(Model *mp, MbRandom *rp, Observations *op, NexusFile *sp, bool dpp, bool adm, int wi, double cp) {

	// initialize some variables indicating what type of restaurant we are
	assumingDpp             = dpp;
	assumingAdmixture       = adm;
	individualsInRestaurant = wi;
	concentrationParm       = cp;

	// set up pointers
	modelPtr                = mp;
	ranPtr                  = rp;
	observationsPtr         = op;
	settingsPtr             = sp;
	
	// initialize the log factorial table
	initializeFactorials();
	
	// allocate the patrons
	initializePatrons();

	// add the tables
	initializeTables();
	
	// initialize table with admixture proportions and counts
	initializeAdmixtureProportions();
	
#	if 0
	cout << "pause" << endl;
	getchar();
	cout << "BEGIN printing restaurant information" << endl;
	print();
	printData();
	cout << "END printing restaurant information" << endl;
	cout << "pause" << endl;
	getchar();
#	endif
}

Restaurant::~Restaurant(void) {

	for (int i=0; i<numPatrons; i++)
		delete patrons[i];
	delete [] patrons;
	delete [] lnFactorial;
	for (vector<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
		delete *p;
}

void Restaurant::addTableToRestaurant(Table *tp) {

	tables.push_back( tp );
}

void Restaurant::deleteTable(Table *tp) {

	bool successfullyDeletedTable = false;
	for (vector<Table *>::iterator t=tables.begin(); t != tables.end(); t++)
		{
		if ( (*t) == tp )
			{
			tables.erase( t );
			delete tp;
			successfullyDeletedTable = true;
			break;
			}
		}
	if (successfullyDeletedTable == false)
		Stmessage::warning("Problem deleting table");
}

int Restaurant::getIndexForTable(Table *tp) {

	for (int i=0; i<getNumTables(); i++)
		{
		if ( getTable(i) == tp )
			return i;
		}
	return -1;
}

void Restaurant::initializeAdmixtureProportions(void) {

	// set up vectors for admixture
	if (assumingAdmixture == true && assumingDpp == false)
		{
		// what type of prior are we placing on the mixing proportions?
		if (settingsPtr->getAdmixtureDist() == "fixed")
			admixtureVarianceParm = settingsPtr->getAdmixtureParm1();
		else if (settingsPtr->getAdmixtureDist() == "exponential")
			admixtureVarianceParm = ranPtr->exponentialRv( settingsPtr->getAdmixtureParm1() );
		else if (settingsPtr->getAdmixtureDist() == "uniform")
			admixtureVarianceParm = ranPtr->uniformRv( settingsPtr->getAdmixtureParm1(), settingsPtr->getAdmixtureParm2() );
		else if (settingsPtr->getAdmixtureDist() == "gamma")
			admixtureVarianceParm = ranPtr->gammaRv( settingsPtr->getAdmixtureParm1(), settingsPtr->getAdmixtureParm2() );
			
		// get the upper and lower limits for the admixture variance parameter
		if (settingsPtr->getAdmixtureDist() == "uniform")
			{
			admixtureLowerLimit = settingsPtr->getAdmixtureParm1();
			admixtureUpperLimit = settingsPtr->getAdmixtureParm2();
			}
		else if (settingsPtr->getAdmixtureDist() == "exponential")
			{
			admixtureLowerLimit = 0.001;
			admixtureUpperLimit = ( 1.0 / settingsPtr->getAdmixtureParm1() ) * 100.0;
			}
		else if (settingsPtr->getAdmixtureDist() == "gamma")
			{
			admixtureLowerLimit = 0.001;
			admixtureUpperLimit = ( settingsPtr->getAdmixtureParm1() / settingsPtr->getAdmixtureParm2() ) * 100.0;
			}
			
		// initialize the counts and probabilities for the mixing proportions
		popCounts = new MbVector<int>[observationsPtr->getNumIndividuals()];
		popProbs  = new MbVector<double>[observationsPtr->getNumIndividuals()];
		for (int i=0; i<observationsPtr->getNumIndividuals(); i++)
			{
			popCounts[i] = MbVector<int>(modelPtr->getNumPopulations());
			popProbs[i] = MbVector<double>(modelPtr->getNumPopulations());
			MbVector<double> alp( modelPtr->getNumPopulations() );
			for (int j=0; j<modelPtr->getNumPopulations(); j++)
				alp[j] = admixtureVarianceParm / modelPtr->getNumPopulations();
			ranPtr->dirichletRv(alp, popProbs[i]);
			for (int j=0; j<modelPtr->getNumPopulations(); j++)
				popCounts[i][j] = 0;
			}
		
		for (int i=0; i<getNumTables(); i++)
			{
			Table *tp = getTable(i);
			for (int j=0; j<tp->getNumPatrons(); j++)
				{
				Patron *p = tp->getPatron(j);
				popCounts[ p->getIndividual() ][ i ]++;
				}
			}
			
		}
}

void Restaurant::initializeFactorials(void) {

	int maxNumFactorialsNeeded = 2 * observationsPtr->getNumIndividuals() * observationsPtr->getNumLoci();
	lnFactorial = new double[maxNumFactorialsNeeded];
	lnFactorial[0] = 0.0;
	for (int n=1; n<maxNumFactorialsNeeded; n++)
		{
		double x = lnFactorial[n-1];
		lnFactorial[n] = x + log( (double)n );
		}
}

void Restaurant::initializePatrons(void) {

	// how many patrons?
	if (assumingAdmixture == false)
		numPatrons = observationsPtr->getNumIndividuals();
	else if (assumingAdmixture == true && individualsInRestaurant == -1)
		numPatrons = observationsPtr->getNumIndividuals() * observationsPtr->getTotalNumAlleleCopies();
	else if (assumingAdmixture == true && individualsInRestaurant >= 0)
		numPatrons = observationsPtr->getTotalNumAlleleCopies();
		
	// allocate the patrons
	patrons = new Patron*[numPatrons];

	// set the index information for each patron
	if (assumingAdmixture == false)
		{
		for (int i=0; i<numPatrons; i++)
			{
			patrons[i] = new PatronIndividual( observationsPtr, i, observationsPtr->getNumLoci(), observationsPtr->getNumUniqueAlleles() );
			patrons[i]->setIndex(i);
			patrons[i]->setIndividual(i);
			}
		}
	else if (assumingAdmixture == true && individualsInRestaurant == -1)
		{
		for (int i=0, m=0; i<observationsPtr->getNumIndividuals(); i++)
			{
			for (int j=0; j<observationsPtr->getNumLoci(); j++)
				{
				int ploidy = 1;
				if (observationsPtr->getIsDiploid(j) == true)
					ploidy = 2;
				for (int k=0; k<ploidy; k++)
					{
					patrons[m] = new PatronAllele( observationsPtr, m, i, j, k, observationsPtr->getNumUniqueAllelesAtLocus(j) );
					m++;
					}
				}
			}
		}
	else if (assumingAdmixture == true && individualsInRestaurant >= 0)
		{
		for (int j=0, m=0; j<observationsPtr->getNumLoci(); j++)
			{
			int ploidy = 1;
			if (observationsPtr->getIsDiploid(j) == true)
				ploidy = 2;
			for (int k=0; k<ploidy; k++)
				{
				patrons[m] = new PatronAllele( observationsPtr, m, individualsInRestaurant, j, k, observationsPtr->getNumUniqueAllelesAtLocus(j) );
				m++;
				}
			}
		}
		
	// set the data for each patron
	for (int i=0; i<numPatrons; i++)
		{
		Patron *p = patrons[i];
		int elemIdx = p->getIndividual();
		p->setData( observationsPtr->getPtrToIndividual(elemIdx) );
		}

	// set the ln likelihood for each data element when placed in a population by itself
	for (int i=0; i<numPatrons; i++)
		{
		Patron *p = patrons[i];
		double x = lnLikelihoodNewTable(p);
		p->setLnLikelihoodAlone(x);
		}		
}

void Restaurant::initializeTables(void) {

	if (assumingDpp == true)
		{
		// Dirichlet process prior model
		for (int i=0; i<numPatrons; i++)
			{
			double probNewTable = concentrationParm / (i + concentrationParm);
			double u = ranPtr->uniformRv();
			bool successfullyAddedPatron = false;
			if (u < probNewTable)
				{
				Table *tab = new Table(this, ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), lnFactorial);
				tables.push_back( tab );
				tab->addPatronToTable( patrons[i] );
				patrons[i]->setTable(tab);
				successfullyAddedPatron = true;
				}
			else
				{
				double sum = probNewTable;
				for (vector<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
					{
					sum += (double)(*p)->getNumPatrons() / (i + concentrationParm);
					if (u < sum)
						{
						(*p)->addPatronToTable( patrons[i] );
						patrons[i]->setTable( (*p) );
						successfullyAddedPatron = true;
						break;
						}
					}
				}
			if (successfullyAddedPatron == false)
				{
				Stmessage::error("Failed to add patron to table");
				return;
				}
			}
		}
	else
		{
		// fixed number of tables
		for (int i=0; i<modelPtr->getNumPopulations(); i++)
			tables.push_back( new Table(this, ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), lnFactorial) );
		MbVector<int> patronIndices(numPatrons);
		for (int i=0; i<numPatrons; i++)
			patronIndices[i] = i;
		for (int i=0; i<tables.size(); i++)
			{
			int whichPatron = (int)(ranPtr->uniformRv()*(numPatrons-i));
			Patron *p = patrons[ patronIndices[whichPatron] ];
			tables[i]->addPatronToTable( p );
			p->setTable(tables[i]);
			int tempIdx = patronIndices[numPatrons-i-1];
			patronIndices[numPatrons-i-1] = patronIndices[whichPatron];
			patronIndices[whichPatron] = tempIdx;
			}
		for (int i=tables.size(); i<numPatrons; i++)
			{
			int whichPatron = (int)(ranPtr->uniformRv()*(numPatrons-i));
			Patron *p = patrons[ patronIndices[whichPatron] ];
			int whichTable = (int)(ranPtr->uniformRv()*tables.size());
			tables[whichTable]->addPatronToTable( p );
			p->setTable(tables[whichTable]);
			int tempIdx = patronIndices[numPatrons-i-1];
			patronIndices[numPatrons-i-1] = patronIndices[whichPatron];
			patronIndices[whichPatron] = tempIdx;
			}
		}
}

double Restaurant::lnLikelihood(void) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnL = 0.0;
	for (vector<Table *>::iterator r=tables.begin(); r != tables.end(); r++)
		lnL += (*r)->lnLikelihood();
	return lnL;
}

double Restaurant::lnLikelihoodNewTable(Patron *p) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnP = 0.0;
	int *x = p->getDataPtr();
	if ( p->getPatronIsIndividual() == true )
		{
		// calculate likelihood for all loci
		for (int l=0; l<observationsPtr->getNumLoci(); l++)
			{
			int k = observationsPtr->getNumUniqueAllelesAtLocus(l);
			int sum1 = 0;
			for (int j=0; j<k; j++)
				{
				lnP += ( lnFactorial[(*x)] );
				sum1 += (*x);
				x++;
				}
			lnP += ( lnFactorial[k - 1] - lnFactorial[k + sum1 - 1] );
			}
		}
	else
		{
		// calculate likelihood for one locus only
		int sum1 = 0;
		for (int j=0; j<p->getNumDataElements(); j++)
			{
			lnP += ( lnFactorial[(*x)] );
			sum1 += (*x);
			x++;
			}
		lnP += ( lnFactorial[p->getNumDataElements() - 1] - lnFactorial[p->getNumDataElements() - sum1 - 1] );
		}
	return lnP;
}

double Restaurant::lnPriorProbability(void) {

	if (assumingAdmixture == false && assumingDpp == false)
		{
		// fixed number of populations with no admixture
		return -numPatrons * log(getNumTables());
		}
	else if (assumingAdmixture == true && assumingDpp == false)
		{
		// fixed number of poulations with admixture
		double lnP = 0.0;
		for (int idx=0; idx<numPatrons; idx++)
			{
			Patron *p = patrons[idx];
			Table *tablePtr = p->getTable();
			lnP += log(popProbs[ p->getIndividual() ][ getIndexForTable(tablePtr) ]);
			}
		return lnP;
		}
	else if (assumingAdmixture == false && assumingDpp == true)
		{
		// number of populations is a random variable, but no admixture
		/* f({\mathbf z}, k | \alpha, c) = \alpha^k {\prod_{i=1}^k(\eta_i - 1)! \over \prod_{i=1}^c(\alpha + i - 1)} */
		double lnSum1 = 0.0;
		for (vector<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
			{
			lnSum1 += lnFactorial[ (*p)->getNumPatrons() - 1 ];
			}
		double lnSum2 = 0.0;
		for (int i=1; i<=numPatrons; i++)
			lnSum2 += log(concentrationParm + i - 1.0);
		double lnP = getNumTables() * log(concentrationParm) + lnSum1 - lnSum2;
		return lnP;
		}
	return 0.0;
}

void Restaurant::normalizeVector(vector<double> &v) {

	// find maximum value
	double lnC = *(max_element( v.begin(), v.end() ));
		
	// rescale to maximum value
	for (int i=0; i<v.size(); i++)
		v[i] -= lnC;
	
	// perform safe exponentation and compute sum for normalizing
	double sum = 0.0;
	for (int i=0; i<v.size(); i++)
		{
		if ( v[i] < -300.0 )
			v[i] = 0.0;
		else
			v[i] = exp( v[i] );
		sum += v[i];
		}
	
	// normalize
	for (int i=0; i<v.size(); i++)
		v[i] /= sum;
}

void Restaurant::print(void) {

	int i = 1;
	for (vector<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
		{
		if ((*p)->getMenuItem() == NULL)
			cout << setw(4) << i++ << ": ";
		else
			cout << setw(4) << i++ << " (" << (*p)->getMenuItem() << "): ";
		(*p)->print();
		}
}

void Restaurant::printData(void) {

	int i = 1;
	for (vector<Table *>::iterator p=tables.begin(); p != tables.end(); p++)
		{
		cout << "Table " << setw(4) << i++ << endl;
		(*p)->printData();
		}
}

int Restaurant::multinomialRv(vector<double> &v) {

	// we assume that v is a probability vector
	int whichCell = 0;
	double u = ranPtr->uniformRv(), sum = 0.0;
	for (vector<double>::iterator prob=v.begin(); prob != v.end(); prob++)
		{
		sum += (*prob);
		if ( u < sum )
			break;
		whichCell++;
		}
	return whichCell;
}

void Restaurant::updateAdmixtureProportions(void) {

	int k = getNumTables();
	MbVector<double> a(k);
	for (int i=0; i<observationsPtr->getNumIndividuals(); i++)
		{
		for (int j=0; j<k; j++)
			a[j] = (popCounts[i][j] + admixtureVarianceParm);
		//cout << "a(" << i << ") = " << fixed << setprecision(10) << a << " : ";
		//cout << popProbs[i] << " -> ";
		ranPtr->dirichletRv( a, popProbs[i] );
		//cout << popProbs[i] << endl;
		}
}

void Restaurant::updateAdmixtureVarianceParm(void) {

	if ( settingsPtr->getAdmixtureDist() == "fixed" )
		return;

	int k = getNumTables();
	MbVector<double> a1(k);
	MbVector<double> a2(k);
	
	// propose a new value for the parameter
	double oldA = admixtureVarianceParm;
	double window = 0.3;
	double newA = oldA + (ranPtr->uniformRv() * window - 0.5 * window);
	if (newA < admixtureLowerLimit)
		newA = 2.0 * admixtureLowerLimit - newA;
	else if (newA > admixtureUpperLimit)
		newA = 2.0 * admixtureUpperLimit - newA;
	
	// calculate the log likelihood ratio
	for (int j=0; j<k; j++)
		{
		a1[j] = oldA / k;
		a2[j] = newA / k;
		}
	double lnR = 0.0;
	for (int i=0; i<observationsPtr->getNumIndividuals(); i++)
		lnR += ( ranPtr->lnDirichletPdf(a2, popProbs[i]) - ranPtr->lnDirichletPdf(a1, popProbs[i]) );
		
	// add in the prior ratio 
	if (settingsPtr->getAdmixtureDist() == "gamma")
		lnR += ( ranPtr->lnGammaPdf(settingsPtr->getAdmixtureParm1(), settingsPtr->getAdmixtureParm2(), newA) - ranPtr->lnGammaPdf(settingsPtr->getAdmixtureParm1(), settingsPtr->getAdmixtureParm2(), oldA) );
	else if (settingsPtr->getAdmixtureDist() == "exponential")
		lnR += ( ranPtr->lnExponentialPdf(settingsPtr->getAdmixtureParm1(), newA) - ranPtr->lnExponentialPdf(settingsPtr->getAdmixtureParm1(), oldA) );
		
	// calculate the acceptance probability
	double r = 0.0;
	if (lnR < -300.0)
		r = 0.0;
	else if (lnR > 0.0)
		r = 1.0;
	else
		r = exp(lnR);
	
	// accept or reject the proposed state
	if (ranPtr->uniformRv() < r)
		admixtureVarianceParm = newA;
	admixtureVarianceParm = 0.05;
}

void Restaurant::updateConcParm(void) {

	updateConcParm( settingsPtr->getConcentrationParm1(), settingsPtr->getConcentrationParm2());
}

void Restaurant::updateConcParm(double a, double b) {

	if (settingsPtr->getConcentrationDist() != "gamma")
		return;

	/* set the parameters */
	int k = getNumTables();
	int n = getNumPatronsInRestaurant();
	double oldAlpha = concentrationParm;
	
	/* Step 1: Draw a Beta(oldAlpha+1, n) distribution to get eta */
	MbVector<double> z(2);
	MbVector<double> f(2);
	z[0] = oldAlpha + 1.0;
	z[1] = (double)n;
	ranPtr->dirichletRv(z, f);
	double eta = f[0];
	
	/* Step 2: Draw a new value for alpha, based on k and eta */
	double u = ranPtr->uniformRv();
	double x = ( a + (double)k - 1.0 ) / ( (double)n * (b - log(eta)) );
	double newAlpha;
	if ( (u / (1.0 - u)) < x)
		newAlpha = ranPtr->gammaRv(a + k, b - log(eta));
	else
		newAlpha = ranPtr->gammaRv(a + k - 1.0, b - log(eta)); // set alpha in the restaurant ...
	concentrationParm = newAlpha;
}

void Restaurant::updateSeatingForFixedTables(void) {

	// allocate memory for reassignment probabilities
	vector<double> probs;
	
	// loop over all patrons, attempting to reseat each
	for (int idx=0; idx<numPatrons; idx++)
		{
		// get a pointer to the element and its current table
		Patron *p = patrons[idx];
		Table *tablePtr = p->getTable();

		// remove the element from the table
		tablePtr->removePatronFromTable(p);
		if (assumingAdmixture == true)
			popCounts[ p->getIndividual() ][ getIndexForTable(tablePtr) ]--;
		
		// calculate the likelihood when the element is seated at all possible tables
		for (int i=0; i<getNumTables(); i++)
			{
			Table *tp = getTable(i);
			double lnC = 0.0;
			if (assumingAdmixture == true)
				lnC = log( popProbs[ p->getIndividual() ][ i ] );
			probs.push_back( tp->lnLikelihood(p) + lnC );
			}
		normalizeVector( probs );

		// pick a table for reassignment
		tablePtr = getTable( multinomialRv(probs) );
			
		// add the element to that table
		tablePtr->addPatronToTable(p);
		if (assumingAdmixture == true)
			popCounts[ p->getIndividual() ][ getIndexForTable(tablePtr) ]++;
			
		// clear the probs vector
		probs.clear();
		}
}

void Restaurant::updateSeatingForDpp(void) {

	for (int idx=0; idx<numPatrons; idx++)
		{
		// get a pointer to the element and its current table
		Patron *p = patrons[idx];
		Table *tablePtr = p->getTable();

		// remove the element from the table
		tablePtr->removePatronFromTable(p);
		if (tablePtr->getNumPatrons() == 0)
			deleteTable(tablePtr);
		
		// calculate the likelihood when the element is seated at all possible tables
		vector<double> probs;
		for (int i=0; i<getNumTables(); i++)
			{
			Table *tp = getTable(i);
			probs.push_back( tp->lnLikelihood(p) + log(tp->getNumPatrons()) );
			}
		probs.push_back( p->getLnLikelihoodalone() + log(concentrationParm) );
		normalizeVector( probs );

		// pick a table for reassignment
		int whichTable = multinomialRv(probs);
			
		// add the element to that table
		if (whichTable < getNumTables())
			{
			tablePtr = getTable(whichTable);
			tablePtr->addPatronToTable(p);
			}
		else
			{
			tablePtr = new Table(this, ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), lnFactorial);
			addTableToRestaurant( tablePtr );
			tablePtr->addPatronToTable(p);
			}
		}
}

#pragma mark Franchise

Franchise::Franchise(Model *mp, MbRandom *rp, Observations *op, NexusFile *sp) {

	// set up pointers
	modelPtr        = mp;
	ranPtr          = rp;
	observationsPtr = op;
	settingsPtr     = sp;
	
	// set the concentration parameters
	restaurantAlpha = modelPtr->getAlpha2();
	menuAlpha       = modelPtr->getAlpha1();
		
	// instantiate restaurants, one for each individual
	for (int i=0; i<observationsPtr->getNumIndividuals(); i++)
		restaurants.push_back( new Restaurant( modelPtr, ranPtr, observationsPtr, settingsPtr, true, true, i, restaurantAlpha ) );
		
	// set up the tables, that represent the global franchise menus
	int tblCnt = 0;
	for (vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		{
		for (int t=0; t<(*p)->getNumTables(); t++)
			{
			Table *restaurantTable = (*p)->getTable(t);
			double probNewTable = menuAlpha / (tblCnt + menuAlpha);
			double u = ranPtr->uniformRv();
			if (u < probNewTable)
				{
				Menu *newMenuItem = new Menu( this, ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), (*p)->getFactorialPtr() );
				newMenuItem->addTableToMenu(restaurantTable);
				menus.push_back( newMenuItem );
				restaurantTable->setMenuItem(newMenuItem);
				}
			else
				{
				double sum = probNewTable;
				for (int i=0; i<menus.size(); i++)
					{
					sum += menus[i]->getNumTables() / (tblCnt + menuAlpha);
					if (u < sum)
						{
						menus[i]->addTableToMenu(restaurantTable);
						restaurantTable->setMenuItem(menus[i]);
						break;
						}
					}
				}
			
			tblCnt++;
			}
		}
		
	// Set up the stirling numbers of the first kind, which are used to calculate the probability
	// of K tables in a restaurant. We are assuming that there is a common concentration parameter
	// for all of the restaurants, which means we cannot perform the simple Gibb's sampling on
	// alpha described by Escobar and West (1995).
	stirlingNums = new Stirling( restaurants[0]->getNumPatronsInRestaurant() );
}

Franchise::~Franchise(void) {

	for (vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		delete *p;
	for (vector<Menu *>::iterator p=menus.begin(); p != menus.end(); p++)
		delete *p;
	delete stirlingNums;
}

void Franchise::deleteMenuItem(Menu *mToDelete) {

	bool successfullyDeletedTable = false;
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		{
		if ( (*m) == mToDelete )
			{
			menus.erase( m );
			delete mToDelete;
			successfullyDeletedTable = true;
			break;
			}
		}
	if (successfullyDeletedTable == false)
		Stmessage::warning("Problem deleting menu table");
}

int Franchise::getIndexForMenuItem(Menu *mp) {

	int menuIdx = 0;
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		{
		if ( (*m) == mp )
			break;
		menuIdx++;
		}
	return menuIdx;
}

int Franchise::getNumPatronsInFranchise(void) {

	int numP = 0;
	for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
		numP += (*r)->getNumPatronsInRestaurant();
	return numP;
}

int Franchise::getNumTablesInFranchise(void) {

	int nt = 0;
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		nt += (*m)->getNumTables();
	return nt;
}

double Franchise::lnLikelihood(void) {

#	if defined(NO_DATA)
	return 0.0;
#	endif

	double lnL = 0.0;
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		lnL += (*m)->lnLikelihood();
	return lnL;
}

double Franchise::lnPriorProbability(void) {
	
	double alpha = restaurantAlpha;
	double gamma = menuAlpha;
		double lnP = 0.0;
	for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
		{
		double lnSum1 = 0.0;
		for (int i=0; i<(*r)->getNumTables(); i++)
			{
			Table *p = (*r)->getTable(i);
			lnSum1 += ranPtr->lnGamma(p->getNumPatrons());
			}
		double lnSum2 = 0.0;
		for (int i=1; i<=(*r)->getNumPatronsInRestaurant(); i++)
			lnSum2 += log(alpha + i - 1.0);
		lnP += (*r)->getNumTables() * log(alpha) + lnSum1 - lnSum2;
		}

	double lnSum1 = 0.0;
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		lnSum1 += ranPtr->lnGamma((*m)->getNumTables());
	double lnSum2 = 0.0;
	for (int i=1; i<=getNumTablesInFranchise(); i++)
		lnSum2 += log(gamma + i - 1.0);
	lnP += getNumMenuItems() * log(gamma) + lnSum1 - lnSum2;
	return lnP;
	
	return 0.0;
}

void Franchise::print(void) {

	int i = 0;
	for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
		{
		cout << "Restaurant " << i++ << ":" << endl;
		(*r)->print();
		}
	cout << "Menu (" << menus.size() << "): ";
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		(*m)->print();
}

void Franchise::printMenuData(void) {

	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		(*m)->printData();
}

int Franchise::numRestaurantTables(void) {

	int nt = 0;
	for (vector<Restaurant *>::iterator p=restaurants.begin(); p != restaurants.end(); p++)
		nt += (*p)->getNumTables();
	return nt;
}

double Franchise::safeExponentiation(double lnR) {

	if (lnR > 0.0)
		return 1.0;
	else if (lnR < -300.0)
		return 0.0;
	else
		return exp(lnR);
}

void Franchise::updateConcParm(void) {

	// update the concentration parameter for the menu items using a Gibb's update
	if (settingsPtr->getConcentrationDist() == "gamma")
		{
		/* set the parameters */
		int k = getNumMenuItems();
		int n = getNumTablesInFranchise();
		double a = settingsPtr->getConcentrationParm1();
		double b = settingsPtr->getConcentrationParm2();
		double oldAlpha = menuAlpha;
		
		/* Step 1: Draw a Beta(oldAlpha+1, n) distribution to get eta */
		MbVector<double> z(2);
		MbVector<double> f(2);
		z[0] = oldAlpha + 1.0;
		z[1] = (double)n;
		ranPtr->dirichletRv(z, f);
		double eta = f[0];
		
		/* Step 2: Draw a new value for alpha, based on k and eta */
		double u = ranPtr->uniformRv();
		double x = ( a + (double)k - 1.0 ) / ( (double)n * (b - log(eta)) );
		double newAlpha;
		if ( (u / (1.0 - u)) < x)
			newAlpha = ranPtr->gammaRv(a + k, b - log(eta));
		else
			newAlpha = ranPtr->gammaRv(a + k - 1.0, b - log(eta)); // set alpha in the restaurant ...
		menuAlpha = newAlpha;
		}
	
	// Update using Metropolis-Hastings sampling. We assume a common concentration parameter
	// for all of the restaurants.
	if (settingsPtr->getAdmixtureConcentrationDist() == "gamma")
		{
		// propose new value
		double oldAlpha = restaurantAlpha;
		double tuning = log(4.0);
		double newAlpha = oldAlpha * exp( tuning*(ranPtr->uniformRv()-0.5) );
		double a = settingsPtr->getAdmixtureConcentrationParm1();
		double b = settingsPtr->getAdmixtureConcentrationParm2();
		
		// get likelihood ratio
		double lnLikelihoodRatio = 0.0;
		stirlingNums->calculateLnProbK(newAlpha);
		for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
			lnLikelihoodRatio += stirlingNums->getLnProbK( (*r)->getNumTables()-1 );
		stirlingNums->calculateLnProbK(oldAlpha);
		for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
			lnLikelihoodRatio -= stirlingNums->getLnProbK( (*r)->getNumTables()-1 );
			
		// get prior ratio
		double lnPriorRatio = ranPtr->lnGammaPdf(a, b, newAlpha) - ranPtr->lnGammaPdf(a, b, oldAlpha);
		
		// get proposal ratio
		double lnProposalRatio = log(newAlpha) - log(oldAlpha);
		
		// acceptance probability
		double r = safeExponentiation(lnLikelihoodRatio + lnPriorRatio + lnProposalRatio);
		
		// accept/reject
		if (ranPtr->uniformRv() < r)
			{
			restaurantAlpha = newAlpha;
			for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
				(*r)->setConcentrationParm(restaurantAlpha);
			}
		}
}

void Franchise::updateSeating(void) {

	// loop over the restaurants
	for (vector<Restaurant *>::iterator r=restaurants.begin(); r != restaurants.end(); r++)
		{
		// loop over the patrons in a restaurant
		for (int idx=0; idx<(*r)->getNumPatronsInRestaurant(); idx++)
			{
			// get a pointer to the patron and its current table
			Patron *p = (*r)->getRestaurantPatron(idx);
			Table *tablePtr = p->getTable();
			Menu *menuItemPtr = tablePtr->getMenuItem();

			// remove the element from the table
			tablePtr->removePatronFromTable(p);
			if (tablePtr->getNumPatrons() == 0)
				{
				menuItemPtr->deleteTableFromMenu(tablePtr);
				if (menuItemPtr->getNumTables() == 0)
					deleteMenuItem(menuItemPtr);
				(*r)->deleteTable(tablePtr);
				}
				
			// calculate the probability of reseating patron at each of the tables
			int n = (*r)->getNumTables() + 1;
			vector<double> probs;
			for (int i=0; i<(*r)->getNumTables(); i++)
				{
				Table *t = (*r)->getTable(i);
				Menu *m = t->getMenuItem();
				int n_jt = t->getNumPatrons();
				probs.push_back( m->lnLikelihood(p) + log(n_jt) );
				}
			vector<double> lnProbsNew;
			probs.push_back( log(restaurantAlpha) + lnLikelihoodNewTable(p, lnProbsNew) );
			(*r)->normalizeVector(probs);
			
			// pick a table
			int whichTable = (*r)->multinomialRv( probs );
			Table *newTablePtr = NULL;
			if  (whichTable < n-1)
				newTablePtr = (*r)->getTable(whichTable);
			
			// seat the patron
			if (newTablePtr == NULL)
				{
				// we are seating the patron at a new table, which may have an already used menu item, or a brand new menu item
				Table *newRestaurantTablePtr = new Table( (*r), ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), (*r)->getFactorialPtr() );
				newRestaurantTablePtr->addPatronToTable( p );
				(*r)->addTableToRestaurant( newRestaurantTablePtr );
				(*r)->normalizeVector(lnProbsNew);
				whichTable = (*r)->multinomialRv( lnProbsNew );
				Menu *newMenuItemPtr = NULL;
				if (whichTable < lnProbsNew.size()-1)
					{
					newMenuItemPtr = menus[whichTable];
					}
				else
					{
					newMenuItemPtr = new Menu( this, ranPtr, observationsPtr, observationsPtr->getNumUniqueAlleles(), (*r)->getFactorialPtr() );
					menus.push_back( newMenuItemPtr );
					}
				newMenuItemPtr->addTableToMenu( newRestaurantTablePtr );
				newRestaurantTablePtr->setMenuItem(newMenuItemPtr);
				}
			else
				{
				// otherwise, we are seating the element at one of the pre-existing tables
				newTablePtr->addPatronToTable( p );
				}
			
		//print();
		//getchar();
			}
		}
}

double Franchise::lnLikelihoodNewTable(Patron *p, vector<double> &pv) {

	// fill in the log probabilities
	vector<double> lnProbs;
	int numTablesInAllRestaurants = numRestaurantTables();
	for (vector<Menu *>::iterator m=menus.begin(); m != menus.end(); m++)
		{
		double x = (double)(*m)->getNumTables() / (numTablesInAllRestaurants + menuAlpha);
		double y = (*m)->lnLikelihood(p);
		double z = log(x) + y;
		lnProbs.push_back( z );
		}
	double z = log( menuAlpha / (numTablesInAllRestaurants + menuAlpha) ) + p->getLnLikelihoodalone();
	lnProbs.push_back( z );
	pv = lnProbs;
	
	// perform the sum, with appropriate scaling
	double lnP = sumLogProbs( lnProbs );
	return lnP;
}

double Franchise::sumLogProbs(vector<double> &v) {

	// find maximum value
	double lnC = *(max_element( v.begin(), v.end() ));
		
	// rescale to maximum value
	for (int i=0; i<v.size(); i++)
		v[i] -= lnC;
	
	// perform safe exponentation and compute sum for normalizing
	double sum = 0.0;
	for (int i=0; i<v.size(); i++)
		{
		double x;
		if ( v[i] < -300.0 )
			x = 0.0;
		else
			x = exp( v[i] );
		sum += x;
		}
	
	// return sum
	return log(sum) + lnC;
}

#pragma mark Model

Model::Model(int idx, NexusFile &settings) {

	/* initialize some variables to default values */
	alpha1         = 1.0;
	conc1Scale     = 1.0;
	conc1Shape     = 1.0;
	conc1Lambda    = 1.0;
	numPopulations = 1;
	numElements    = 1;

	/* set the model index */
	modelIndex = idx; 
	
	/* get a pointer to the observations object */
	observationsPtr = settings.samplePtr();
	
	/* set up the random number object */
	ranPtr = new MbRandom( (long int)(time(NULL) + modelIndex) );

	/* determine which model we have from the following potential models:
	
	  Model  Number of populations  Admixture
	  ---------------------------------------
	  M1     Fixed                  No 
	  M2     Fixed                  Yes
	  M3     Random Variable        No
	  M4     Random Variable        Yes
	  --------------------------------------- */            
	modelId = -1;
	if ( settings.getNumPops() > 0 && settings.getAdmixture() == false )
		modelId = 1;
	else if ( settings.getNumPops() > 0 && settings.getAdmixture() == true)
		modelId = 2;
	else if ( settings.getNumPops() == 0 && settings.getAdmixture() == false )
		modelId = 3;
	else if ( settings.getNumPops() == 0 && settings.getAdmixture() == true )
		modelId = 4;
	if (modelId == -1)
		{
		Stmessage::warning("Unknown model settings");
		return;
		}
		
	/* some variables that are only relevant for reporting the model that is being used */
	string numPopulationsModel = "";
	concentrationParm1Model = "";
	concentrationParm2Model = "";
	admixtureModel = "";
	string admixtureYnModel = "No";
	char tempStr[100];
		
	/* set the relevant parameters depending on the model 
	
	   1. set the number of populations (if the number of populations is fixed) */
	if (modelId == 1 || modelId == 2)
		{
		numPopulations = settings.getNumPops();
		sprintf(tempStr, "Fixed(%d)", settings.getNumPops());
		numPopulationsModel = tempStr;
		}
	else
		{
		sprintf(tempStr, "Random Variable (Dirichlet Process Prior)");
		numPopulationsModel = tempStr;
		}
		
	/* set the number of elements to arrange among tables */
	if (modelId == 1 || modelId == 3)
		{
		numElements = settings.getNumIndividuals();
		numElementsPerRestaurant = settings.getNumIndividuals();
		}
	else if (modelId == 2)
		{
		numElements = settings.getNumIndividuals() * observationsPtr->getTotalNumAlleleCopies();
		numElementsPerRestaurant = settings.getNumIndividuals() * observationsPtr->getTotalNumAlleleCopies();
		}
	else if (modelId == 4)
		{
		numElements = settings.getNumIndividuals() * observationsPtr->getTotalNumAlleleCopies();
		numElementsPerRestaurant = observationsPtr->getTotalNumAlleleCopies();
		}
		
	/* 2. set the concentration parameter controlling the number of populations */
	if (modelId == 3 || modelId == 4)                                                              
		{
		if (settings.getConcentrationDist() == "fixed")
			{
			alpha1 = settings.getConcentrationParm1();
			sprintf(tempStr, "Fixed(%1.3lf)", alpha1);
			concentrationParm1Model = tempStr;
			}
		else if (settings.getConcentrationDist() == "fixed_expk")
			{
			alpha1 = calcAlphaFromExpectedNumberOfTables( numElements, settings.getConcentrationParm1() );
			sprintf(tempStr, "Fixed(E(K)=%1.3lf)", settings.getConcentrationParm1());
			concentrationParm1Model = tempStr;
			}
		else if (settings.getConcentrationDist() == "exponential")
			{
			conc1Lambda = settings.getConcentrationParm1();
			alpha1 = ranPtr->exponentialRv( conc1Lambda );
			sprintf(tempStr, "Exponential(%1.3lf)", conc1Lambda);
			concentrationParm1Model = tempStr;
			}
		else if (settings.getConcentrationDist() == "gamma")
			{
			conc1Shape = settings.getConcentrationParm1();
			conc1Scale = settings.getConcentrationParm2();
			alpha1 = ranPtr->gammaRv( conc1Shape, conc1Scale );
			sprintf(tempStr, "Gamma(%1.3lf,%1.3lf)", conc1Shape, conc1Scale);
			concentrationParm1Model = tempStr;
			}
		}
		
	/* 3. set the concentration parameter controlling the number of mixture components for each individual */
	if (modelId == 4)                                                                              	
		{
		admixtureYnModel = "Yes";
		if (settings.getAdmixtureConcentrationDist() == "fixed")
			{
			alpha2 = settings.getAdmixtureConcentrationParm1();
			sprintf(tempStr, "Fixed(%1.3lf)", alpha2);
			concentrationParm2Model = tempStr;
			}
		else if (settings.getAdmixtureConcentrationDist() == "fixed_expk")
			{
			alpha2 = calcAlphaFromExpectedNumberOfTables( numElementsPerRestaurant*observationsPtr->getNumIndividuals(), settings.getAdmixtureConcentrationParm1() );
			sprintf(tempStr, "Fixed(E(K)=%1.3lf)", settings.getAdmixtureConcentrationParm1());
			concentrationParm2Model = tempStr;
			}
		else if (settings.getAdmixtureConcentrationDist() == "exponential")
			{
			conc2Lambda = settings.getAdmixtureConcentrationParm1();
			alpha2 = ranPtr->exponentialRv( conc2Lambda );
			sprintf(tempStr, "Exponential(%1.3lf)", conc2Lambda);
			concentrationParm2Model = tempStr;
			}
		else if (settings.getAdmixtureConcentrationDist() == "gamma")
			{
			conc2Shape = settings.getAdmixtureConcentrationParm1();
			conc2Scale = settings.getAdmixtureConcentrationParm2();
			alpha2 = ranPtr->gammaRv( conc2Shape, conc2Scale );
			sprintf(tempStr, "Gamma(%1.3lf,%1.3lf)", conc2Shape, conc2Scale);
			concentrationParm2Model = tempStr;
			}
		}
		
	/* 4. set the parameter controlling admixture proportions */
	if (modelId == 2)
		{
		admixtureYnModel = "Yes";
		if (settings.getAdmixtureDist() == "fixed")
			{
			alpha3 = settings.getAdmixtureParm1();
			sprintf(tempStr, "Fixed(%1.3lf)", alpha3);
			admixtureModel = tempStr;
			}
		else if (settings.getAdmixtureDist() == "exponential")
			{
			admLambda = settings.getAdmixtureParm1();
			alpha3 = ranPtr->exponentialRv(admLambda);
			sprintf(tempStr, "Exponential(%1.3lf)", admLambda);
			admixtureModel = tempStr;
			}
		else if (settings.getAdmixtureDist() == "gamma")
			{
			admShape = settings.getAdmixtureParm1();
			admScale = settings.getAdmixtureParm2();
			alpha3 = ranPtr->gammaRv( admShape, admScale );
			sprintf(tempStr, "Gamma(%1.3lf,%1.3lf)", admShape, admScale);
			admixtureModel = tempStr;
			}
		else if (settings.getAdmixtureDist() == "uniform")
			{
			admUniform[0] = settings.getAdmixtureParm1();
			admUniform[1] = settings.getAdmixtureParm2();
			alpha3 = ranPtr->uniformRv(admUniform[0], admUniform[1]);
			sprintf(tempStr, "Uniform(%1.3lf,%1.3lf)", admUniform[0], admUniform[1]);
			admixtureModel = tempStr;
			}
		}
	
	/* print information to the screen about the model settings */
	if (modelIndex == 0)
		{
		stprint("   Model\n");
		if (modelId == 1)
			{
			stprint("      Number of populations  = %s\n", numPopulationsModel.c_str() );
			stprint("      Admixture              = %s\n", admixtureYnModel.c_str() );
			}
		else if (modelId == 2)
			{
			stprint("      Number of populations        = %s\n", numPopulationsModel.c_str() );
			stprint("      Admixture                    = %s\n", admixtureYnModel.c_str() );
			stprint("      Admixture variance parameter = %s\n", admixtureModel.c_str() );
			}
		else if (modelId == 3)
			{
			stprint("      Number of populations       = %s\n", numPopulationsModel.c_str() );
			stprint("      DPP concentration parameter = %s\n", concentrationParm1Model.c_str() );
			stprint("      Admixture                   = %s\n", admixtureYnModel.c_str() );
			}
		else if (modelId == 4)
			{
			stprint("      Number of populations             = %s\n", numPopulationsModel.c_str() );
			stprint("      DPP concentration parameter       = %s\n", concentrationParm1Model.c_str() );
			stprint("      Admixture                         = %s\n", admixtureYnModel.c_str() );
			stprint("      Admixture concentration parameter = %s\n", concentrationParm2Model.c_str() );
			}
		}

	/* set up the restaurant or franchise, that holds the tables */
	restaurant = NULL, franchise = NULL;
	if (modelId == 1)
		restaurant = new Restaurant(this, ranPtr, observationsPtr, &settings, false, false, -1, 0.0);
	else if (modelId == 2)
		restaurant = new Restaurant(this, ranPtr, observationsPtr, &settings, false, true, -1, 0.0);
	else if (modelId == 3)
		restaurant = new Restaurant(this, ranPtr, observationsPtr, &settings, true, false, -1, getAlpha1());
	else if (modelId == 4)
		franchise = new Franchise(this, ranPtr, observationsPtr, &settings);
}

Model::~Model(void) {

	delete ranPtr;
	if (restaurant != NULL)
		delete restaurant;
	if (franchise != NULL)
		delete franchise;
}

double Model::calcAlphaFromExpectedNumberOfTables(int n, double expT) {

	double a = 0.000001;
	double ea = expNumTables(n, a);
	bool goUp;
	if (ea < expT)
		goUp = true;
	else
		goUp = false;
	double increment = 0.1;
	while ( fabs(ea - expT) > 0.000001 )
		{
		if (ea < expT && goUp == true)
			{
			a += increment;
			}
		else if (ea > expT && goUp == false)
			{
			a -= increment;
			}
		else if (ea < expT && goUp == false)
			{
			increment /= 2.0;
			goUp = true;
			a += increment;
			}
		else
			{
			increment /= 2.0;
			goUp = false;
			a -= increment;
			}
		ea = expNumTables(n, a);
		}
	//stprint(ea << " <-> " << expT << " " << "alpha=" << a \n;
	return a;
}

double Model::expNumTables(int n, double a) {

	double expectedNum = 0.0;
	for (int i=1; i<=n; i++)
		expectedNum += ( 1.0 / (i - 1.0 + a) );
	expectedNum *= a;
	return expectedNum;
}

double Model::expNumTables(int n, double alpha, double beta) {

	MbVector<double>cats(5000);
	ranPtr->discretizeGamma(cats, alpha, beta, cats.size(), false);
	
	double ave = 0.0;
	for (int k=0; k<cats.size(); k++)
		{
		double a = cats[k];
		double expectedNum = 0.0;
		for (int i=1; i<=n; i++)
			expectedNum += ( 1.0 / (i - 1.0 + a) );
		expectedNum *= a;
		ave += expectedNum;
		}
	ave /= cats.size();
	
	return ave;
}

bool Model::isNumPopsRandom(void) {

	if ( modelId == 1 || modelId == 2 )
		return false;
	return true;
}

double Model::getAdmixtureVarianceParm(void) {

	return restaurant->getAdmixtureVarianceParm();
}

int Model::getCurrentNumberOfPopulations(void) {

	if (modelId == 1 || modelId == 2 || modelId == 3)
		return restaurant->getNumTables();
	else if (modelId == 4)
		return franchise->getNumMenuItems();
	return 0;
}

double Model::lnLikelihood(void) {

	if (modelId == 1 || modelId == 2 || modelId == 3)
		return restaurant->lnLikelihood();
	else if (modelId == 4)
		return franchise->lnLikelihood();
	return 0.0;
}

double Model::lnPriorProbability(void) {

	if (modelId == 1 || modelId == 2 || modelId == 3)
		return restaurant->lnPriorProbability();
	else if (modelId == 4)
		return franchise->lnPriorProbability();
	return 0.0;
}

PopPartition* Model::getPartition(void) {

	PopPartition *partPtr;
	if (modelId == 4)
		partPtr = new PopPartition( franchise );
	else
		partPtr = new PopPartition( restaurant );
	return partPtr;
}

void Model::print(void) {

	if (modelId == 4)
		franchise->print();
	else
		restaurant->print();
}

void Model::updateSeating(void) {

	if (modelId == 1 || modelId == 2)
		restaurant->updateSeatingForFixedTables();
	else if (modelId == 3)
		restaurant->updateSeatingForDpp();
	else if (modelId == 4)
		franchise->updateSeating();
}

void Model::updateAdmixtureProportions(void) {

	if (modelId == 2)
		restaurant->updateAdmixtureProportions();
}

void Model::updateAdmixtureVarianceParm(void) {

	if (modelId == 2)
		restaurant->updateAdmixtureVarianceParm();
}

void Model::updateConcParm(void) {

	if (modelId == 3)
		restaurant->updateConcParm();
	else if (modelId == 4)
		franchise->updateConcParm();
}

