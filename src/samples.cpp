#include "iomanager.h"
#include "model.h"
#include "nexusfile.h"
#include "observation.h"
#include "partition.h"
#include "samples.h"
#include "st.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>

#if defined(MAC_GUI)
#import <Cocoa/Cocoa.h>
#import "AppController.h"
#endif




Sample::Sample(void) {

	admixtureVariance       = -1.0;
	admixtureVariancePrior  = "";
	cycle                   = -1;
	dppAlphaPopulation      = -1.0;
	dppAlphaPopulationPrior = "";
	dppAlphaAdmixture       = -1.0;
	dppAlphaAdmixturePrior  = "";
	partition               = NULL;
}

Sample::~Sample(void) {

	if (partition != NULL)
		delete partition;
}

string Sample::getSampleString(void) {

	string s = "";
	char tempC[50];
	
	// add cycle 
	sprintf(tempC, "%d\t", cycle);
	s += tempC;
	
	// add likelihood
	sprintf(tempC, "%1.4f\t", lnLike);
	s+= tempC;
	
	// add variance parameter
	if (admixtureVariance >= 0.0)
		{
		sprintf(tempC, "%1.4lf\t", admixtureVariance);
		s += tempC;
		}
	else if (dppAlphaPopulation >= 0.0)
		{
		sprintf(tempC, "%1.4lf\t", dppAlphaPopulation);
		s += tempC;
		if (dppAlphaAdmixture >= 0.0)
			{
			sprintf(tempC, "%1.4lf\t", dppAlphaAdmixture);
			s += tempC;
			}
		}
		
	// add partition
	if (partition != NULL)
		{
		for (int i=0; i<partition->getNumElements(); i++)
			{
			sprintf(tempC, "%d\t", partition->getElement(i));
			s += tempC;
			}
		}
	
	return s;
}

McmcSamples::McmcSamples(void) {

	meanPartition               = NULL;
	jointProbs                  = NULL;
	numSamplesMeanPartIsBasedOn = 0;
	bestMeanPartScore           = 0.0;
	gpsInfoPresent              = false;
	meanPartitionMixtureProbs   = NULL;
}

McmcSamples::~McmcSamples(void) {

	deleteSamples();
}

void McmcSamples::calcMarginalLike(NexusFile &settings) {

	// check that we have at least one sample
	if (isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// put together a list of likelihoods
	vector<double> likes;
	for (int i=0; i<getNumSamples(); i++)
		{
		if (i >= settings.getBurnInMarginalLike())
			likes.push_back( samples[i]->getLnLike() );
		}

	// check that we have at least one partition in the list
	if (likes.size() == 0)
		{
		Stmessage::warning("All of the samples have been discarded, using the specified burn-in");
		return;
		}
		
	// calculate the marginal likelihood
	int n = 0;
	double scaler = 0.0;
	double a = 0.0;
	double harmonicMean;
	for (vector<double>::iterator p=likes.begin(); p != likes.end(); p++)
		{
		double y = (*p);
		y = -y;
		if (n == 0)
			{
			scaler = y;
			}
		else
			{
			if (y > scaler)
				{
				double diff = scaler - y;
				if (diff < -300.0)
					a = 0.0;
				else
					a *= exp(diff);
				scaler = y;
				}
			}
		y -= scaler;
		double z = exp(y);
		a += z;
		n++;
		harmonicMean = log( (double)n ) - (scaler + log(a));
		}
		
	stprint("   Marginal likelihood calculated using the method of Newton & Raftery (1994),\n");
	stprint("   based on %d MCMC samples:\n\n", likes.size() );
	stprint("   Marginal log likelihood = %1.2lf\n\n", harmonicMean );
	
	stprint("   Newton, M. A., and A. E. Raftery. 1994. Approximate Bayesian inference\n");
	stprint("      by the weighted likelihood bootstrap. Journal of the Royal Statistical\n");
	stprint("      Society, Series B. 3:3-48.\n");
}

void McmcSamples::calcMeanPart(NexusFile &settings) {

	// check that we have at least one sample
	if (isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// put together a list of partitions
	vector<PopPartition *> partList;
	for (int i=0; i<getNumSamples(); i++)
		{
		if (i >= settings.getBurnInMeanPart())
			partList.push_back( samples[i]->getPartition() );
		}
		
	// check that we have at least one partition in the list
	if (partList.size() == 0)
		{
		Stmessage::warning("All of the samples have been discarded, using the specified burn-in");
		return;
		}

	// open file
	ofstream fileOut;
	IoManager *io = settings.getMeanPartOutPtr();
	if (settings.getSaveToFileMeanPart() == true)
		{
		if ( io->testDirectory() == false )
			{
			Stmessage::error("Cannot find directory for output of mean partition");
			return;
			}
		if ( io->testFile() == true && settings.getWarnUser() == true)
			{
			if ( MyString::queryUser("   A file with the name \"" + io->getFileName() + "\" already exists.\n   Do you with to overwrite the file?") == false )
				{
				stprint("   Terminating Showmeanpartition\n");
				return;
				}
			}
		if ( io->openFile(fileOut) == false )
			{
			Stmessage::error("Could not open file for output of mean partition");
			return;
			}
		}
		
	// check that we haven't already made the mean partition
	if ( meanPartition != NULL && numSamplesMeanPartIsBasedOn != partList.size() )
		{
		delete meanPartition;
		meanPartition = NULL;
		}
		
#	if defined (MAC_GUI)
	[appControllerPtr convertStatusWindowIndeterminate];
	[appControllerPtr turnOffAllToolBarItems];
	[appControllerPtr setProcessRunning:YES];
#	endif

	// make the mean partition
	bool cancelMeanPartitionProcess = false;
	if (meanPartition == NULL)
		meanPartition = new PopPartition(partList, bestMeanPartScore, var, cancelMeanPartitionProcess);
		
	// check to see if the mean partition process was canceled
	if ( cancelMeanPartitionProcess == true )
		{
		numSamplesMeanPartIsBasedOn = 0;
		delete meanPartition;
		meanPartition = NULL;
#		if defined(MAC_GUI)
		stprint("   Canceled heuristic search for mean partition\n");
		[appControllerPtr resetStatusWindow];
		[appControllerPtr updateToolBar];
		[appControllerPtr setProcessRunning:NO];
#		endif
		return;
		}
		
	// we can set the sample size of the current mean partition
	numSamplesMeanPartIsBasedOn = partList.size();
	
	// display the mean partition
	if (settings.getSaveToFileMeanPart() == true)
		fileOut << "Identifier" << '\t' << "Population" << endl;
	stprint("   Sum of squares score for mean partition = %1.2lf\n", bestMeanPartScore );
	stprint("   Mean partition based on %d MCMC samples", numSamplesMeanPartIsBasedOn );
	stprint("   Mean partition:\n\n");
	stprint("   %8s    %4s     %4s\n", "Identifier", "Population", "Variance");
	stprint("   -------------------------------------\n");
	for (int i=0; i<meanPartition->getNumElements(); i++)
		{
		PatronId *pid = patronIds[i];
		char tempC[50];
		if (pid->isIndividual == true)
			sprintf(tempC, "%d", pid->individualId);
		else
			sprintf(tempC, "(%d,%d,%d)", pid->individualId, pid->locusId, pid->alleleId);
		string tempStr = tempC;
		stprint("   %10s -- %10d   %10.3lf\n", tempStr.c_str(), meanPartition->getElement(i), var[i] );
		if (settings.getSaveToFileMeanPart() == true)
			fileOut << tempStr << '\t' << meanPartition->getElement(i) << endl;
		}
	stprint("   -------------------------------------\n");
	
	// close file
	if (settings.getSaveToFileTogetherness() == true)
		io->closeFile(fileOut);
				
	// set up a vector with the mixture probabilities
	int n = getNumPatrons();
	int d = getDegreeOfMeanPartition();
	meanPartitionMixtureProbs = new float[n * d];
	for (int i=0; i<n*d; i++)
		meanPartitionMixtureProbs[i] = 0.0;
	float *f = &meanPartitionMixtureProbs[0];
	for (int i=1; i<=n; i++)
		{
		getMixtureProbsForMeanPartition( i, f, d );
		f += d;
		}
	if (n != meanPartition->getNumElements())
		{
		// only print if we have admixture, in which case we show the admixture proportions
		stprint("\n");
		stprint("   Mean Partition for individuals, showing admixture proportions:\n\n");
		stprint("   %10s", "Identifier");
		for (int j=0; j<d; j++)
			stprint("%6d", j+1);
		stprint("\n");
		stprint("   --------------");
		for (int j=0; j<d; j++)
			stprint("------");
		stprint("\n");
		f = &meanPartitionMixtureProbs[0];
		for (int i=1; i<=n; i++)
			{
			stprint("   %10d -- ", i);
			for (int j=0; j<d; j++)
				stprint("%6.3f", f[j]);
			stprint("\n");
			f += d;
			}
		stprint("   --------------");
		for (int j=0; j<d; j++)
			stprint("------");
		stprint("\n");
		}
		
#	if defined(MAC_GUI)
	[appControllerPtr resetStatusWindow];
	[appControllerPtr updateToolBar];
	[appControllerPtr setProcessRunning:NO];
#	endif
}

void McmcSamples::calcTogetherness(NexusFile &settings) {

	// check that we have at least one sample
	if (isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// put together a list of partitions
	vector<PopPartition *> partList;
	for (int i=0; i<getNumSamples(); i++)
		{
		if (i >= settings.getBurnInTogetherness())
			partList.push_back( samples[i]->getPartition() );
		}
		
	// check that we have at least one partition in the list
	if (partList.size() == 0)
		{
		Stmessage::warning("All of the samples have been discarded, using the specified burn-in");
		return;
		}
		
	// count the number of elements per locus
	int numAllelesPerIndividual = 0;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( (*p)->individualId == 1 )
			numAllelesPerIndividual++;
		}
		
	// count the number of individuals
	int numIndividuals = patronIds.size() / numAllelesPerIndividual;

	// open file
	ofstream fileOut;
	IoManager *io = settings.getTogethernessOutPtr();
	if (settings.getSaveToFileTogetherness() == true)
		{
		if ( io->testDirectory() == false )
			{
			Stmessage::error("Cannot find directory for output of togetherness probabilities");
			return;
			}
		if ( io->testFile() == true && settings.getWarnUser() == true )
			{
			if ( MyString::queryUser("   A file with the name \"" + io->getFileName() + "\" already exists.\n   Do you with to overwrite the file?") == false )
				{
				stprint("   Terminating Showtogetherness\n");
				return;
				}
			}
		if ( io->openFile(fileOut) == false )
			{
			Stmessage::error("Could not open file for output of togetherness probabilities");
			return;
			}
		}

	// calculate togetherness probabilities
#	if defined (MAC_GUI)
	[appControllerPtr convertStatusWindowIndeterminate];
#	endif
	if (jointProbs != NULL)
		{
		delete [] jointProbs;
		jointProbs = NULL;
		}
	if (jointProbs == NULL)
		jointProbs = new float[numIndividuals * (numIndividuals-1) / 2];
	if (settings.getSaveToFileTogetherness() == true)
		fileOut << "First Individual" << '\t' << "Second Individual" << '\t' << "Pr[Individuals in same population]" << endl;
	stprint("   Probability that pairs of individuals are from the same\n");
	stprint("   population, based on %d MCMC samples:\n\n", partList.size() );
	for (int i=0; i<3; i++)
		stprint("   %10s%7s", "Pair", "Prob." );
	stprint("\n");
	for (int i=0; i<3; i++)
		stprint("   -----------------");
	stprint("\n");
	int cnt = 0;
	for (int i=0, k=0; i<numIndividuals; i++)
		{
		int posI = i * numAllelesPerIndividual;
		for (int j=i+1; j<numIndividuals; j++)
			{
			int posJ = j * numAllelesPerIndividual;
			double prob = 0.0;
			for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
				{
				int a = posI, b = posJ, numSame = 0;
				for (int k=0; k<numAllelesPerIndividual; k++)
					{
					if ( (*p)->getElement(a) == (*p)->getElement(b) )
						numSame++;
					}
				prob += (double)numSame / numAllelesPerIndividual;
				}
			char tempC[50];
			sprintf(tempC, "(%d,%d)", i+1, j+1);
			string tempStr = tempC;
			stprint("   %10s%7.3lf", tempStr.c_str(), prob / partList.size() );
			if (settings.getSaveToFileTogetherness() == true)
				fileOut << i+1 << '\t' << j+1 << '\t' << fixed << setprecision(3) << prob / partList.size() << endl;
			cnt++;
			if (cnt % 3 == 0)
				stprint("\n");
			jointProbs[k++] = prob / partList.size();
			}
		}
	if (cnt % 3 != 0)
		stprint("\n");
	for (int i=0; i<3; i++)
		stprint("   -----------------");
	stprint("\n");
#	if defined (MAC_GUI)
	[appControllerPtr resetStatusWindow];
#	endif

	// close file
	if (settings.getSaveToFileTogetherness() == true)
		io->closeFile(fileOut);
}

void McmcSamples::calcNumPops(NexusFile &settings) {

	// check that we have at least one sample
	if (isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// put together a list of partitions
	vector<PopPartition *> partList;
	for (int i=0; i<getNumSamples(); i++)
		{
		if (i >= settings.getBurnInNumPops())
			partList.push_back( samples[i]->getPartition() );
		}
		
	// check that we have at least one partition in the list
	if (partList.size() == 0)
		{
		Stmessage::warning("All of the samples have been discarded, using the specified burn-in");
		return;
		}

	// open file
	ofstream fileOut;
	IoManager *io = settings.getNumPopsOutPtr();
	if (settings.getSaveToFileNumPops() == true)
		{
		if ( io->testDirectory() == false )
			{
			Stmessage::error("Cannot find directory for output of number of populations distribution");
			return;
			}
		if ( io->testFile() == true && settings.getWarnUser() == true )
			{
			if ( MyString::queryUser("   A file with the name \"" + io->getFileName() + "\" already exists.\n   Do you with to overwrite the file?") == false )
				{
				stprint("   Terminating Shownumpops\n");
				return;
				}
			}
		if ( io->openFile(fileOut) == false )
			{
			Stmessage::error("Could not open file for output of number of populations distribution");
			return;
			}
		}

	// find the maximum degree
	int maxDegree = 0;
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		if ( (*p)->getDegree() > maxDegree )
			maxDegree = (*p)->getDegree();
		}
	
	// count the number of samples for each degree
	int *cnts = new int[maxDegree];
	for (int i=0; i<maxDegree; i++)
		cnts[i] = 0;
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		cnts[ (*p)->getDegree()-1 ]++;
		
	// display the results
	if (settings.getSaveToFileNumPops() == true)
		fileOut << "Num. Pops." << '\t' << "Prob." << endl;
	stprint("   Probability distribution for the number of populations,\n");
	stprint("   based on %d MCMC samples:\n\n", partList.size() );
	stprint("   %10s    %6s\n", "Num. Pops.", "Prob." );
	stprint("   --------------------\n");
	for (int i=0; i<maxDegree; i++)
		{
		stprint("   %10d -- %6.2lf\n", i+1, (double)cnts[i] / partList.size() );
		if (settings.getSaveToFileNumPops() == true)
			fileOut << i+1 << '\t' << fixed << setprecision(2) << (double)cnts[i] / partList.size() << endl;
		}
	stprint("   --------------------\n");
	
	// close file
	if (settings.getSaveToFileNumPops() == true)
		io->closeFile(fileOut);
}

void McmcSamples::deleteSamples(void) {

	if (meanPartition != NULL)
		delete meanPartition;
	meanPartition = NULL;
	if (meanPartitionMixtureProbs != NULL)
		delete [] meanPartitionMixtureProbs;
	meanPartitionMixtureProbs = NULL;
	if (jointProbs != NULL)
		delete [] jointProbs;
	for (vector<Sample *>::iterator s=samples.begin(); s != samples.end(); s++)
		delete (*s);
	samples.clear();
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		delete (*p);
	patronIds.clear();
	header.clear();
	numSamplesMeanPartIsBasedOn = 0;
	bestMeanPartScore = 0.0;
	gpsInfoPresent = false;
}

string McmcSamples::getHeaderString(Model *mp) {

	string s = "";

	// add cycle 
	s += "Cycle\t";
	
	// add likelihood
	s+= "lnL\t";
	
	// add variance parameter
	if (mp->getModelId() == 2)
		s += "v~" + mp->getAdmixtureModel() + '\t';
	else if (mp->getModelId() == 3 || mp->getModelId() == 4)
		{
		s += "a1~" + mp->getConcentrationParm1Model() + '\t';
		if (mp->getModelId() == 4)
			s += "a2~" + mp->getConcentrationParm2Model() + '\t';
		}
		
	// add the header information
	gpsInfoPresent = false;
	if (mp->getModelId() == 4)
		{
		Franchise *fp = mp->getFranchisePtr();
		for (int r=0; r<fp->getNumRestaurants(); r++)
			{
			Restaurant *rp = fp->getRestaurant(r);
			for (int i=0; i<rp->getNumPatronsInRestaurant(); i++)
				{
				Patron *pp = rp->getRestaurantPatron(i);
				char tempC[50];
				sprintf(tempC, "I(%d,%d,%d", pp->getIndividual()+1, pp->getLocus()+1, pp->getAllele()+1);
				s += tempC;
				if (pp->hasGpsCoordinates() == true)
					{
					gpsInfoPresent = true;
					sprintf(tempC, ",{%1.6lf,%1.6lf}", pp->getLatitude(), pp->getLongitude() );
					s += tempC;
					}
				s += ")\t";
				}
			}
		}
	else
		{
		Restaurant *rp = mp->getRestaurantPtr();
		for (int i=0; i<rp->getNumPatronsInRestaurant(); i++)
			{
			Patron *pp = rp->getRestaurantPatron(i);
			char tempC[50];
			if (pp->getPatronIsIndividual() == false)
				sprintf(tempC, "I(%d,%d,%d", pp->getIndividual()+1, pp->getLocus()+1, pp->getAllele()+1);
			else
				sprintf(tempC, "I(%d", pp->getIndividual()+1);
				s += tempC;
			if (pp->hasGpsCoordinates() == true)
				{
				gpsInfoPresent = true;
				sprintf(tempC, ",{%1.6lf,%1.6lf}", pp->getLatitude(), pp->getLongitude() );
				s += tempC;
				}
			s += ")\t";
			}
		}
	
	return s;
}

double McmcSamples::getLatitude(int i) {

	int n = 0;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( hasGps( (*p) ) == true )
			{
			if ( (*p)->isIndividual == true )
				{
				if (n == i)
					return (*p)->latitude;
				n++;
				}
			else
				{
				if ( (*p)->locusId == 0 && (*p)->alleleId == 0 )
					{
					if (n == i)
						return (*p)->latitude;
					n++;
					}
				}
			}
		}
	return 0.0;
}

double McmcSamples::getLongitude(int i) {

	int n = 0;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( hasGps( (*p) ) == true )
			{
			if ( (*p)->isIndividual == true )
				{
				if (n == i)
					return (*p)->longitude;
				n++;
				}
			else
				{
				if ( (*p)->locusId == 0 && (*p)->alleleId == 0 )
					{
					if (n == i)
						return (*p)->longitude;
					n++;
					}
				}
			}
		}
	return 0.0;
}

int McmcSamples::getNumGpsCoordinates(void) {

	int n = 0;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( hasGps( (*p) ) == true )
			{
			if ( (*p)->isIndividual == true )
				{
				n++;
				}
			else
				{
				if ( (*p)->locusId == 0 && (*p)->alleleId == 0 )
					n++;
				}
			}
		}
	return n;
}

int McmcSamples::getNumPatrons(void) {

	vector<int> individualsFound;
	for(vector<PatronId *>::iterator p=patronIds.begin(); p!=patronIds.end(); p++)
		{
		int indIdx = (*p)->individualId;
		vector<int>::iterator location = find( individualsFound.begin(), individualsFound.end(), indIdx );
		if ( location == individualsFound.end() )
			individualsFound.push_back( indIdx );
		}
	return individualsFound.size();
}

int McmcSamples::getNumMixturesForPatronOnMeanPartition(int indIdx) {

	vector<int> mixturesFound;
	int i = 0;
	for(vector<PatronId *>::iterator p=patronIds.begin(); p!=patronIds.end(); p++)
		{
		int x = (*p)->individualId;
		if ( x == indIdx )
			{
			int popIdx = meanPartition->getElement(i);
			vector<int>::iterator location = find( mixturesFound.begin(), mixturesFound.end(), popIdx );
			if ( location == mixturesFound.end() )
				mixturesFound.push_back( popIdx );
			}
		i++;
		}
	return mixturesFound.size();
}

int McmcSamples::getSmallestPatronId(void) {

	int smallestPatronId = 1000000000;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( (*p)->individualId < smallestPatronId )
			smallestPatronId = (*p)->individualId;
		}
	return smallestPatronId;
}

int McmcSamples::getLargestPatronId(void) {

	int largestPatronId = -1000000000;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		if ( (*p)->individualId > largestPatronId )
			largestPatronId = (*p)->individualId;
		}
	return largestPatronId;
}

int McmcSamples::getDegreeOfMeanPartition(void) {

	return meanPartition->getDegree();
}

float* McmcSamples::getMeanPartitionMixtureProbs(void) {

	return meanPartitionMixtureProbs;
}

void McmcSamples::getMixtureProbsForMeanPartition(int indIdx, float *probs, int numMixtures) {

	for (int j=0; j<numMixtures; j++)
		probs[j] = 0.0;
		
	int i = 0;
	for (vector<PatronId *>::iterator p=patronIds.begin(); p != patronIds.end(); p++)
		{
		int x = (*p)->individualId;
		if ( x == indIdx )
			{
			int popIdx = meanPartition->getElement(i) - 1;
			if (popIdx < 0 || popIdx >= numMixtures)
				cout << "ERROR: Problem getting mixture probabilities" << endl;
			probs[popIdx] += 1.0;
			}
		i++;
		}

	float sum = 0.0;
	for (int j=0; j<numMixtures; j++)
		sum += probs[j];
	for (int j=0; j<numMixtures; j++)
		probs[j] /= sum;
}

bool McmcSamples::hasGps(PatronId *p) {

	if ( p->latitude > -90.0 && p->latitude < 90.0 && p->longitude > -180.0 && p->longitude < 180.0 )
		return true;
	return false;
}

int McmcSamples::getNumIndividuals(void) {

	return samples[0]->getPartition()->getNumElements();
}

void McmcSamples::interpretInd(string h, int &a, int &b, int &c, bool &d, double &e, double &f) {

	// set all of the variables to a base state
	a = 0;
	b = 0;
	c = 0;
	d = false;
	e = 1000.0;
	f = 1000.0;
	
	// check that the word contains a complete description
	int cnt = 0;
	for (int i=0; i<h.size(); i++)
		{
		if (h[i] == '(' || h[i] == '{')
			cnt++;
		else if (h[i] == ')' || h[i] == '}')
			cnt--;
		}
	if (cnt != 0)
		{
		cerr << "   ERROR: Unable to interpret header information" << endl;
		exit(1);
		}
		
	// read the information in the description
	string s = "";
	cnt = 0;
	bool readingGps = false;
	for (int i=0; i<h.size(); i++)
		{
		if ( h[i] == ',' || h[i] == ')' || h[i] == '}' )
			{
			if (s == "")
				continue;
			//cout << s << endl;
			istringstream buf(s);
			if (readingGps == false)
				{
				int v;
				buf >> v;
				if (cnt == 0)
					a = v;
				else if (cnt == 1)
					b = v;
				else if (cnt == 2)
					c = v;
				else
					{
					cerr << "   ERROR: Unable to interpret indiviudal information" << endl;
					exit(1);
					}
				cnt++;
				}
			else
				{
				double v;
				buf >> v;
				if (h[i] == ',')
					e = v;
				else
					f = v;
				}
			s = "";
			}
		else if ( h[i] != 'i' && h[i] != '(' && h[i] != '{' && h[i] != ' ' )
			{
			s += h[i];
			}
		if ( h[i] == '{' )
			readingGps = true;
		else if ( h[i] == '}' )
			readingGps = false;
		}
	if (cnt > 1)
		d = false;
	else 
		d = true;
}

bool McmcSamples::isEmpty(void) {

	if (samples.size() == 0)
		return true;
	return false;
}

bool McmcSamples::readFile(IoManager &mngr, NexusFile &settings) {

	// perhaps clear out any samples already in memory
	if (isEmpty() == false)
		{
		if ( settings.getWarnUser() == true )
			if ( MyString::queryUser("   Overwrite MCMC samples currently in memory?") == false )
				return false;
		stprint("   Deleting samples from memory\n");
		deleteSamples();
		}

	// open the file
	ifstream fileIn;
	if ( mngr.openFile(fileIn) == false )
		{
		Stmessage::error("Could not open input file");
		return false;
		}
		
	// read the file
	string linestring = "";
	int line = 0;
	stprint("   Reading samples from file %s\n", mngr.getFileName().c_str());
	gpsInfoPresent = false;
	while( getline(fileIn, linestring).good() )
		{
		istringstream linestream(linestring);
		int ch;
		string word = "";
		int column = 0;
		Sample *newSample = NULL;
		vector<int> rgf;
		do
			{
			word = "";
			linestream >> word;
			MyString::lowerCase(word);
			
			if (word != "")
				{
				if (line == 0)
					{
					// read the header
					addHeader( word );
					if (word[0] == 'i' && word[1] == '(')
						{
						int a, b, c;
						bool d;
						double e, f;
						interpretInd( word, a, b, c, d, e, f );
						PatronId *pid = new PatronId();
						pid->individualId = a;
						pid->locusId      = b;
						pid->alleleId     = c;
						pid->isIndividual = d;
						pid->latitude     = e;
						pid->longitude    = f;
						if ( e >= -90.0 && e <= 90.0 && f >= -180.0 && f <= 180.0 )
							gpsInfoPresent = true;
						addPatronId(pid);
						}
					}
				else
					{
					// read in a sample
					if (newSample == NULL)
						newSample = new Sample();
					string h = header[column];
					istringstream buf(word);
					double v;
					buf >> v;

					if ( h == "cycle" )
						newSample->addCycle( (int)v );
					else if ( h == "lnl" )
						newSample->addLnLike( (float)v );
					else if ( h[0] == 'v' && h[1] == '~' )
						newSample->addAdmixtureVariance( (float)v );
					else if ( h[0] == 'a' && h[1] == '1' && h[2] == '~' )
						newSample->addDppAlphaPopulation( (float)v );
					else if ( h[0] == 'a' && h[1] == '2' && h[2] == '~' )
						newSample->addDppAlphaAdmixture( (float)v );
					else if ( h[0] == 'i' && h[1] == '(' )
						rgf.push_back( (int)v );
					
					}
				//cout << word << " " << column << endl;
				column++;		
				}
				
			} while ( (ch=linestream.get()) != EOF );
		
		if (newSample != NULL)
			{
			PopPartition *newPartition = NULL;
			if (rgf.size() > 0)
				{
				newPartition = new PopPartition(rgf);
				newSample->addPartition(newPartition);
				}
			addSample(newSample);
			newSample = NULL;
			}
		line++;
		}	

	// close the file
	mngr.closeFile(fileIn);
	stprint("      Number of samples read into memory = %d\n", getNumSamples() );
	
	return true;
}

void McmcSamples::setPatronIds(Model *mp) {

	if (mp->getModelId() == 4)
		{
		Franchise *fp = mp->getFranchisePtr();
		for (int r=0; r<fp->getNumRestaurants(); r++)
			{
			Restaurant *rp = fp->getRestaurant(r);
			for (int i=0; i<rp->getNumPatronsInRestaurant(); i++)
				{
				Patron *pp        = rp->getRestaurantPatron(i);
				PatronId *pid     = new PatronId();
				pid->individualId = pp->getIndividual()+1;
				pid->locusId      = pp->getLocus()+1;
				pid->alleleId     = pp->getAllele()+1;
				pid->isIndividual = false;
				pid->latitude     = pp->getLatitude();
				pid->longitude    = pp->getLongitude();
				addPatronId(pid);
				}
			}
		}
	else
		{
		Restaurant *rp = mp->getRestaurantPtr();
		for (int i=0; i<rp->getNumPatronsInRestaurant(); i++)
			{
			Patron *pp = rp->getRestaurantPatron(i);
			PatronId *pid = new PatronId();
			if (pp->getPatronIsIndividual() == false)
				{
				pid->individualId = pp->getIndividual()+1;
				pid->locusId      = pp->getLocus()+1;
				pid->alleleId     = pp->getAllele()+1;
				pid->isIndividual = false;
				pid->latitude     = pp->getLatitude();
				pid->longitude    = pp->getLongitude();
				}
			else
				{
				pid->individualId = pp->getIndividual()+1;
				pid->locusId      = 0;
				pid->alleleId     = 0;
				pid->isIndividual = true;
				pid->latitude     = pp->getLatitude();
				pid->longitude    = pp->getLongitude();
				}
			addPatronId(pid);
			}
		}
}





