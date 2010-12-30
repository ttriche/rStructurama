#include "MbBitfield.h"
#include "observation.h"
#include "st.h"
#include <iostream>
#include <iomanip>
#include <list>
#include <algorithm>

#define MISSING 1010101

using namespace std;



Observations::Observations(int ni, int nl) {

	numIndividuals = ni;
	numLoci = nl;

	for (int i=0; i<numIndividuals; i++)
		individuals.push_back( new Individual(this, numLoci) );
		
	numAllelesAtLocus = new int[numLoci];
	isDiploid = new bool[numLoci];
}

Observations::~Observations(void) {

	for (vector<Individual *>::iterator p=individuals.begin(); p != individuals.end(); p++)
		delete *p;
	delete [] numAllelesAtLocus;
	delete [] isDiploid;
}

void Observations::print(void) {

	int longestName = 11;
	for (int i=0; i<numIndividuals; i++)
		if (individuals[i]->getName().length() > longestName)
			longestName = individuals[i]->getName().length();

	stprint("\n");
	stprint("      %*s -- ", longestName+2, "Individuals" );
	for (int i=0; i<numLoci; i++)
		stprint("%9d ", i+1 );
	stprint("\n");
	stprint("      ---");
	for (int i=0; i<longestName+2; i++)
		stprint("-");
	for (int i=0; i<numLoci; i++)
		stprint("----------");
	stprint("\n");

	for (int i=0; i<numIndividuals; i++)
		{
		Individual *indPtr = individuals[i];
		string word = indPtr->getName();
		stprint("      %*s -- ", longestName+2, word.c_str() );
		for (int j=0; j<numLoci; j++)
			{
			stprint("(");
			if ( indPtr->getNumAllelesAtLocus(j) == 1 )
				{
				int a0 = indPtr->getRawInfoAt(j, 0);
				if (a0 != MISSING)
					stprint("%4d   ", a0 );
				else
					stprint("%4s   ", "?" );
				}
			else
				{
				int a0 = indPtr->getRawInfoAt(j, 0);
				int a1 = indPtr->getRawInfoAt(j, 1);
				if (a0 != MISSING)
					stprint("%3d", a0 );
				else
					stprint("%3s", "?" );
				stprint(",");
				if (a1 != MISSING)
					stprint("%3d", a1 );
				else
					stprint("%3s", "?" );
				}
			stprint(") ");
			}
		if (indPtr->hasGpsCoordinates() == true)
			stprint(" {%1.4lf,%1.4lf}", indPtr->getLatitude(), indPtr->getLongitude());
		stprint("\n");
		}

	stprint("      ---");
	for (int i=0; i<longestName+2; i++)
		stprint("-");
	for (int i=0; i<numLoci; i++)
		stprint("----------");
	stprint("\n");
}

void Observations::printRecoded(void) {

	int longestName = 11;
	for (int i=0; i<numIndividuals; i++)
		if (individuals[i]->getName().length() > longestName)
			longestName = individuals[i]->getName().length();

	stprint("\n");
	stprint("      %*s -- ", longestName+2, "Individuals" );
	for (int i=0; i<numLoci; i++)
		stprint("%9d ", i+1 );
	stprint("\n");
	stprint("      ---");
	for (int i=0; i<longestName+2; i++)
		stprint("-");
	for (int i=0; i<numLoci; i++)
		stprint("----------");
	stprint("\n");

	for (int i=0; i<numIndividuals; i++)
		{
		Individual *indPtr = individuals[i];
		string word = indPtr->getName();
		stprint("      %*s -- ", longestName+2, word.c_str() );
		for (int j=0; j<numLoci; j++)
			{
			stprint("(");
			if ( indPtr->getNumAllelesAtLocus(j) == 1 )
				{
				int a0 = indPtr->getInfoAt(j, 0);
				if (a0 != MISSING)
					stprint("%4d   ", a0 );
				else
					stprint("%4s   ", "?" );
				}
			else
				{
				int a0 = indPtr->getInfoAt(j, 0);
				int a1 = indPtr->getInfoAt(j, 1);
				if (a0 != MISSING)
					stprint("%3d", a0 );
				else
					stprint("%3s", "?" );
				stprint(",");
				if (a1 != MISSING)
					stprint("%3d", a1 );
				else
					stprint("%3s", "?" );
				}
			stprint(") ");
			}
		stprint("\n");
		}

	stprint("      ---");
	for (int i=0; i<longestName+2; i++)
		stprint("-");
	for (int i=0; i<numLoci; i++)
		stprint("----------");
	stprint("\n");
}

list<int> Observations::getListOfAllelesAtLocus(int l) {

	/* get a list of the unique alleles at this locus and check
	   that all of the individuals have the same number of alleles
	   observed at each locus (i.e., that they are all haploid or
	   all diploid */
	list<int> observedAlleles;
	int na = 0;
	for (int i=0; i<numIndividuals; i++)
		{
		Individual *indPtr = individuals[i];
		for (int a=0; a<indPtr->getNumAllelesAtLocus(l); a++)
			{
			int x = indPtr->getRawInfoAt(l, a);
			if (x != MISSING)
				observedAlleles.push_back( x );
			}
		
		if (i == 0)
			na = indPtr->getNumAllelesAtLocus(l);
		else
			{
			if (indPtr->getNumAllelesAtLocus(l) != na)
				{
				stprint("Different number of alleles observed at this locus\n");
				exit(1);
				}
			}
		}
	if (na == 1)
		isDiploid[l] = false;
	else if (na == 2)
		isDiploid[l] = true;
	else
		{
		stprint("Locus is neither haploid nor diploid");
		exit(1);
		}
	observedAlleles.sort();
	observedAlleles.unique();
	return observedAlleles;
}

void Observations::recodeInfo(void) {
	
	/* initialize information for each locus */
	for (int l=0; l<numLoci; l++)
		{
		/* get a list of the unique alleles at this locus */
		list<int> observedAlleles = getListOfAllelesAtLocus(l);
		
		/* set the number of unique alleles observed at the locus */
		numAllelesAtLocus[l] = observedAlleles.size();
		if (numAllelesAtLocus[l] == 0)
			{
			stprint("No unique (non-ambiguous) alleles observed\n" );
			exit(1);
			}
		
		/* print information on the number of unique alleles observed 
		   and also on the ploidy of each locus */
		string ploidy = "";
		if (isDiploid[l] == true)
			ploidy = "diploid";
		else
			ploidy = "haploid";
		stprint( "      Number of unique alleles at locus %d = %d (%s)\n", l+1, numAllelesAtLocus[l], ploidy.c_str() );
		}
		
	/* count up the number of alleles, in total, across all loci */
	numUniqueAlleles = 0;
	for (int l=0; l<numLoci; l++)
		numUniqueAlleles += numAllelesAtLocus[l];
		
	/* count up the number of alleles for each individual */
	totalNumAlleleCopies = 0;
	for (int l=0; l<numLoci; l++)
		{
		totalNumAlleleCopies++;
		if (isDiploid[l] == true)
			totalNumAlleleCopies++;
		}
		
	/* remember where the first bit is for each locus */
	positionForLocus = new int[numLoci];
	int bitPos = 0;
	for (int l=0; l<numLoci; l++)
		{
		positionForLocus[l] = bitPos;
		bitPos += numAllelesAtLocus[l];
		}
	
	/* allocate the information for each individual */
	for (int i=0; i<numIndividuals; i++)
		individuals[i]->allocInfo(numUniqueAlleles);
		
	for (int l=0; l<numLoci; l++)
		{
		list<int> observedAlleles = getListOfAllelesAtLocus(l);
		/* now, initialize the recoded allele information for each individual */
		for (int i=0; i<numIndividuals; i++)
			{
			Individual *indPtr = individuals[i];
			for (int a=0; a<indPtr->getNumAllelesAtLocus(l); a++)
				{
				/* get the original denotation of the allele provided by the user */
				int x = indPtr->getRawInfoAt(l, a);
				
				if ( x == MISSING )
					{
					indPtr->setInfo(l, a, -1);
					continue;
					}
				
				/* find the original denotation of the allele in the 
				   list of the unique alleles at the locus */
				int loc = 0;
				for (list<int>::iterator p=observedAlleles.begin(); p != observedAlleles.end(); p++)
					{
					if ( (*p) == x )
						break;
					loc++;
					}		
				
				/* recode the allele */
				if ( loc == numAllelesAtLocus[l] )
					{
					stprint( "Could not find allele during recoding\n");
					exit(1);
					}
				else
					{
					indPtr->setInfo(l, a, loc);
					}
				}
			/*for (list<int>::iterator p=observedAlleles.begin(); p != observedAlleles.end(); p++)
				stprint( (*p) << " ";
			stprint( \n;*/
			}
		}

#	if 0
	for (int i=0; i<numIndividuals; i++)
		{
		stprint( "individual " << i << ":" << \n;
		individuals[i]->printRecoded();
		}
	printRecoded();
	getchar();
#	endif
}

Individual::Individual(Observations *op, int nl) {

	obsPtr          = op;
	numLoci         = nl;
	isInfoAllocated = false;
	latitude        = 1000.0;
	longitude       = 1000.0;

	rawInfo = new int*[nl];
	rawInfo[0] = new int[nl * 2];
	for (int i=1; i<nl; i++)
		rawInfo[i] = rawInfo[i-1] + 2;
	for (int i=0; i<nl; i++)
		for (int j=0; j<2; j++)
			rawInfo[i][j] = 0;
 	numAlleles = new int[nl];
	for (int i=0; i<nl; i++)
		numAlleles[i] = 0;
}

Individual::~Individual(void) {

	delete [] rawInfo[0];
	delete [] rawInfo;
	delete [] numAlleles;
	if (isInfoAllocated == true)
		freeInfo();
}

void Individual::allocInfo(int x) {

	if (isInfoAllocated == true)
		freeInfo();
		
	for (int i=0; i<2; i++)
		info[i] = new MbBitfield(x);
		
	isInfoAllocated = true;
}

void Individual::freeInfo(void) {

	for (int i=0; i<2; i++)
		delete info[i];
	isInfoAllocated = false;
}

bool Individual::hasGpsCoordinates(void) {

	if (latitude > -90.0 && latitude < 90.0 && longitude > -180.0 && longitude < 180.0)
		return true;
	return false;
}

void Individual::printRecoded(void) {

	for (int a=0; a<2; a++)
		{
		stprint("%d: ", a );
		for (int l=0, k=0; l<numLoci; l++)
			{
			for (int j=0; j<obsPtr->getNumUniqueAllelesAtLocus(l); j++)
				{
				if ( info[a]->isBitSet(k++) )
					stprint("1");
				else
					stprint("0");
				}
			stprint(" ");
			}
		stprint("\n");
		}
}

void Individual::setInfo(int l, int a, int s) {

	if (s < 0)
		return;
	int whichBitToSet = obsPtr->getPositionForLocus(l) + s;
	info[a]->setBit( whichBitToSet );
}

int Individual::getInfoAt(int l, int a) { 

	int whichBitToSet = obsPtr->getPositionForLocus(l);
	int alleleId = -1;
	for (int j=0; j<obsPtr->getNumUniqueAllelesAtLocus(l); j++)
		{
		if ( info[a]->isBitSet(whichBitToSet) == true )
			{
			alleleId = j;
			break;
			}
		whichBitToSet++;
		}
		
	return alleleId;
}

int Individual::getNumUniqueAllelesAtLocus(int l) {

	return obsPtr->getNumUniqueAllelesAtLocus(l);
}

int Individual::getNumUniqueAlleles(void) {

	return obsPtr->getNumUniqueAlleles();
}


