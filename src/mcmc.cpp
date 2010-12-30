#include "iomanager.h"
#include "MbRandom.h"
#include "mcmc.h"
#include "model.h"
#include "nexusfile.h"
#include "partition.h"
#include "samples.h"
#include "st.h"
#include <iostream>
#include <ctime>

#ifdef MAC_GUI
#import <Cocoa/Cocoa.h>
#import "AppController.h"
#endif

using namespace std;

Mcmc::Mcmc(NexusFile &settings) {

	// set the parameters of the Markov chain
	numCycles       = settings.getNgen();
	printFrequency  = settings.getPrintfreq();
	sampleFrequency = settings.getSamplefreq();
	numChains       = settings.getNumChains();
	temperature     = settings.getTemperature();
	mcmcSamplesPtr  = settings.getMcmcSamplesPtr();
	ranPtr          = new MbRandom();
	settingsPtr     = &settings;
	
	// instantiate the model(s)
	models = new Model*[numChains];
	for (int chain=0; chain<numChains; chain++)
		models[chain] = new Model(chain, settings);
}

Mcmc::~Mcmc(void) {

	for (int chain=0; chain<numChains; chain++)
		delete models[chain];
	delete [] models;
	delete ranPtr;
}

void Mcmc::attemptSwap(int **si) {

	/* pick two chains at random */
	int chain1 = (int)(ranPtr->uniformRv() * numChains);
	int chain2 = chain1;
	do
		{
		chain2 = (int)(ranPtr->uniformRv() * numChains);
		} while (chain1 == chain2);
	if (chain1 > chain2)
		{
		int temp = chain1;
		chain1 = chain2;
		chain2 = temp;
		}
		
	/* get information for acceptance prob */
	double lnLike_1 = models[chain1]->lnLikelihood();
	double lnLike_2 = models[chain2]->lnLikelihood();
	double lnPrior_1 = models[chain1]->lnPriorProbability();
	double lnPrior_2 = models[chain2]->lnPriorProbability();
	int index_1 = models[chain1]->getModelIndex();
	int index_2 = models[chain2]->getModelIndex();
	
	/* calculate acceptance probability */
	double lnR = (getHeat(index_2, temperature) * (lnLike_1 + lnPrior_1) + getHeat(index_1, temperature) * (lnLike_2 + lnPrior_2) ) - 
	             (getHeat(index_1, temperature) * (lnLike_1 + lnPrior_1) + getHeat(index_2, temperature) * (lnLike_2 + lnPrior_2) );
	double r = 0.0;
	if (lnR > 0.0)
		r = 1.0;
	else if (lnR < -100.0)
		r = 0.0;
	else
		r = exp(lnR);
		
	/* accept/reject proposed swap */
	bool isAccepted = false;
	if (ranPtr->uniformRv() < r)
		isAccepted = true;
		
	/* update state of chain */
	if ( isAccepted == true )
		{
		models[chain1]->setModelIndex(index_2);
		models[chain2]->setModelIndex(index_1);
		}

	/* update count information */
	int i = index_1;
	int j = index_2;
	if ( index_1 > index_2 )
		{
		i = index_2;
		j = index_1;
		}
	if ( isAccepted == true )
		si[i][j]++;
	si[j][i]++;		
}

int Mcmc::findColdChain(void) {

	int coldId = 0;
	for (int i=0; i<numChains; i++)
		{
		if ( models[i]->getModelIndex() == 0 )
			{
			coldId = i;
			break;
			}
		}
	return coldId;	
}

void Mcmc::runChain(void) {
	
	stprint("   Running chain\n");

	// set up file for output
	string op = settingsPtr->getChainOutPtr()->getFilePath();
	string of = settingsPtr->getChainOutPtr()->getFileName();
	string outputFileAndPath = op + "/" + of;
	bool filePresent;
	outputFile.open( outputFileAndPath.c_str(), ios::out );
	if ( !outputFile )
		filePresent = false;
	else
		filePresent = true;
		
	// warn user if file already exists (the warning occurs on selecting the
	// output file for the Mac GUI version)
#	if !defined(MAC_GUI)
	if (filePresent == true && settingsPtr->getWarnUser() == true)
		{
		if ( MyString::queryUser("   Overwrite output file named \"" + of + "\"?") == false )
			{
			stprint("   Terminating MCMC analysis\n");
			outputFile.close();
			return;
			}
		}
#	endif

	// delete any samples that may be in memory
	if (mcmcSamplesPtr->isEmpty() == false)
		{
		if (settingsPtr->getWarnUser() == true)
			{
			if ( MyString::queryUser("   Replace MCMC samples currently in memory?") == false )
				{
				stprint("   Terminating MCMC analysis\n");
				return;
				}
			}
		mcmcSamplesPtr->deleteSamples();
		stprint("      Deleting samples from memory\n");
		}
		
	// set up a matrix of acceptance rates for the Metropolis coupling
	int **swaps = NULL;
	if (numChains > 1)
		{
		swaps = new int*[numChains];
		swaps[0] = new int[numChains * numChains];
		for (int i=1; i<numChains; i++)
			swaps[i] = swaps[i-1] + numChains;
		for (int i=0; i<numChains; i++)
			for (int j=0; j<numChains; j++)
				swaps[i][j] = 0;
		}
	
#	if defined(MAC_GUI)
	// change the status window to its "I'm doing something" look
	[appControllerPtr convertStatusWindow];
	[appControllerPtr turnOffAllToolBarItems];
	[appControllerPtr setProcessRunning:YES];
	int progressInterval = numCycles / 100;
	if (progressInterval < 1)
		progressInterval = 1;
#	endif

	// run the chain
	startCpuTime = clock();
	bool canceledChain = false;
	for (int n=1; n<=numCycles; n++)
		{
		// cycle through cold & heated chains
		for (int chain=0; chain<numChains; chain++)
			{
			// propose a new seating arrangement (Gibbs)
			models[chain]->updateSeating();
			
			// propose a new concentration parameter (Gibbs)
			models[chain]->updateConcParm();

			// propose new admixture proportions for each individual (Gibbs)
			models[chain]->updateAdmixtureProportions();
			
			// propose a new admixture parameter (variance parameter of Dirichlet, Metropolis-Hastings)
			models[chain]->updateAdmixtureVarianceParm();
			}
			
		// attempt to swap the states of the chain
		if (numChains > 1)
			attemptSwap(swaps);
			
#		if !defined (MAC_GUI)
		// print to the screen
		printToScreen(n);
#		endif

		// sample the chain
		sampleChain(n);
		
#		if defined (MAC_GUI)
		// update the status
		if ( n % progressInterval == 0 || n == numCycles )
			{
			double indicatorProgress = ((double)n / numCycles) * 100.0;
			[appControllerPtr updateProgressBar:indicatorProgress];
			}
		if ( [appControllerPtr cancelAction] == YES )
			{
			canceledChain = true;
			break;
			}
#		endif
		}
	endCpuTime = clock();
	
	// summarize chain statistics
	stprint("   Markov chain completed, taking %1.2f seconds of CPU time\n", (endCpuTime - startCpuTime)/(float)CLOCKS_PER_SEC );

	// print swap summary
	printSwapSummary(swaps);
	
	// close the output file
	outputFile.close();

	if (canceledChain == false)
		stprint("   Successfully ran Markov chain Monte Carlo analysis\n");
	else
		stprint("   Canceled Markov chain Monte Carlo analysis\n");

	// free memory
	if (numChains > 1)
		{
		delete [] swaps[0];
		delete [] swaps;
		}

#	if defined(MAC_GUI)
	// reset the status window to its normal look
	[appControllerPtr resetStatusWindow];
	[appControllerPtr updateToolBar];
	[appControllerPtr setProcessRunning:NO];
#	endif
}

void Mcmc::sampleChain(int n) {

	if (n == 1)
		{
		// add the header information
		string headerStr = mcmcSamplesPtr->getHeaderString( models[findColdChain()] );
		outputFile << headerStr << endl;
		mcmcSamplesPtr->setPatronIds( models[findColdChain()] );
		}

	if (n % sampleFrequency == 0)
		{
		// save the sample
		
		// find the cold chain
		int coldId = findColdChain();
		
		// make a new sample
		Sample *sample = new Sample();
		
		// add cycle information
		sample->addCycle(n);
		
		// add likelihood information
		sample->addLnLike( models[coldId]->lnLikelihood() );
		
		// add the partition
		PopPartition *newPart = models[coldId]->getPartition();
		sample->addPartition(newPart);
		
		// add the variance parameter or the population alpha
		if (models[coldId]->getModelId() == 2)
			{
			sample->addAdmixtureVariance( models[coldId]->getRestaurantPtr()->getAdmixtureVarianceParm() );
			}
		else if (models[coldId]->getModelId() == 3)
			{
			sample->addDppAlphaPopulation( models[coldId]->getRestaurantPtr()->getConcentrationParm() );
			}
		else if (models[coldId]->getModelId() == 4)
			{
			sample->addDppAlphaPopulation( models[coldId]->getFranchisePtr()->getMenuAlpha() );
			sample->addDppAlphaAdmixture( models[coldId]->getFranchisePtr()->getRestaurantAlpha() );
			}
			
		// finally, add the sample to the set of samples
		mcmcSamplesPtr->addSample( sample );
		
		// and print the sample to the file
		string sampleStr = sample->getSampleString();
		outputFile << sampleStr << endl;
		}
}

void Mcmc::printToScreen(int n) {

	if ( n % printFrequency == 0 )
		{
		int coldId = findColdChain();
		clock_t curCpuTime = clock();
		float numSecsElapsed = (curCpuTime - startCpuTime) / (float)CLOCKS_PER_SEC;
		float numSecsToGo    = (numCycles - n) * ((float)numSecsElapsed / n);
		string timeToGo = MyString::formatTime(numSecsToGo);
		stprint("      Generation %d (%s remaining)", n, timeToGo.c_str() );
		if ( numChains > 1 )
			{
			stprint(" [");
			for (int k=0; k<numChains; k++)
				{
				if ( k == coldId )
					stprint("C");
				else
					stprint("*");
				if ( k != numChains - 1 )
					stprint(" ");
				}
			stprint("]");
			}
		stprint("\n");
		}
}

void Mcmc::printSwapSummary(int **si) {

	if ( numChains > 1 )
		{
		int largestSwapEntry = 0;
		for (int i=0; i<numChains; i++)
			for (int j=0; j<numChains; j++)
				if (si[i][j] > largestSwapEntry)
					largestSwapEntry = si[i][j];
		int nDigits = (int)((log((double)largestSwapEntry) / log(10.0))) + 1;
		if (nDigits < 4)
			nDigits = 4;
							
		stprint("\n");
		stprint("   Information on swaps:\n\n");
		stprint("           ");
		for (int i=0; i<numChains; i++)
			stprint("%d ", i+1 );
		stprint("\n");
		stprint("   -------");
		for (int i=0; i<numChains; i++)
			{
			for (int k=0; k<nDigits; k++)
				stprint("-");
			stprint("-");
			}
		stprint("\n");
		for (int i=0; i<numChains; i++)
			{
			stprint("   %4d -- ", i+1 );
			for (int j=0; j<numChains; j++)
				{
				if (i == j)
					{
					for (int k=0; k<nDigits; k++)
						stprint("-");
					stprint(" ");
					}
				else
					{
					if (i < j)
						stprint("%1.2lf ", (double)si[i][j]/si[j][i] );
					else
						stprint("%d ", si[i][j] );
					}
				}
			stprint("\n");
			}
		stprint("\n   Upper diagonal: Proportion of accepted swaps\n");
		stprint("   Lower diagonal: Number of attempted swaps\n\n");
		}
}




