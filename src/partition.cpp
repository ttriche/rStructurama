#include "hungarian.h"
#include "machineInfo.h"
#include "MbBitfield.h"
#include "model.h"
#include "partition.h"
#include "st.h"
#include <ctime>
#include <vector>

#if defined(MAC_GUI)
#import <Cocoa/Cocoa.h>
#import "AppController.h"
#endif

// this might be easier using OpenMP...
#ifdef THREADED_MEAN_PART
#include <pthread.h>

typedef struct {
	int                         elem;
	PopPartition*               partPtr;
	std::vector<PopPartition *> *list;
	double                      curBest;
	int*                        elemIndexPtr;
	double*                     distPtr;
} ThreadSettings;

void *threadDist(void *mstPtr);
#endif

using namespace std;

// Chinese Restaurant Process: find E[rp]
PopPartition::PopPartition(Restaurant *rp) {

	// how many elements in PopPartition?
	numElements = rp->getNumPatronsInRestaurant();
	
	// allocate the PopPartition
	allocatePartition();
	
	// set up the PopPartition
	for (int t=0; t<rp->getNumTables(); t++)
		{
		Table *tp = rp->getTable(t);
		for (int p=0; p<tp->getNumPatrons(); p++)
			{
			Patron *pp = tp->getPatron(p);
			int elementId = pp->getIndex() + 0;
			part[elementId] = t;
			}
		}
	
	// put PopPartition into RGF form
	setRgf();
	
	// get the degree
	setDegree();
	
	// set the subsets
	setSubsets();
}

PopPartition::PopPartition(Franchise *fp) {

	// how many elements in partition?
	numElements = fp->getNumPatronsInFranchise();
	
	// allocate the partition
	allocatePartition();

	// set up the partition
	int offSet = 0;
	for (int r=0; r<fp->getNumRestaurants(); r++)
		{
		Restaurant *rp = fp->getRestaurant(r);
		for (int t=0; t<rp->getNumTables(); t++)
			{
			Table *tp = rp->getTable(t);
			Menu *mp = tp->getMenuItem();
			int menuIdx = fp->getIndexForMenuItem(mp);
			for (int p=0; p<tp->getNumPatrons(); p++)
				{
				Patron *pp = tp->getPatron(p);
				int elementId = pp->getIndex() + offSet;
				part[elementId] = menuIdx;
				}
			}
		offSet += rp->getNumPatronsInRestaurant();
		}

	// put partition into RGF form
	setRgf();

	// get the degree
	setDegree();

	// set the subsets
	setSubsets();
}

PopPartition::PopPartition(PopPartition &p) {

	// how many elements in partition?
	numElements = p.numElements;

	// allocate the partition
	allocatePartition();

	// set up the partition (we assume the copied partition is in RGF form)
	for (int i=0; i<numElements; i++)
		part[i] = p.part[i];
		
	// set the degree
	degree = p.degree;

	// set the subsets
	setSubsets();
}

PopPartition::PopPartition(vector<int> &rgf) {

	// how many elements in partition?
	numElements = rgf.size();

	// allocate the partition
	allocatePartition();

	// set up the partition (we assume the copied partition is in RGF form)
	for (int i=0; i<numElements; i++)
		part[i] = rgf[i];

	// put partition into RGF form
	setRgf();

	// get the degree
	setDegree();

	// set the subsets
	setSubsets();
}

PopPartition::PopPartition(vector<PopPartition *> &partList, double &score, vector<double> &var, bool &cancelMeanPartitionProcess) {

#	ifndef THREADED_MEAN_PART

	// BEGIN serial calculation of mean partition

	// check that all of the partitions in the list have the same size
	int partitionSize = partList[0]->getNumElements();
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		if ( (*p)->getNumElements() != partitionSize )
			{
			cerr << "ERROR: PopPartition list has an inconsistent number of elements" << endl;
			exit(1);
			}
		}
		
	// set the number of elements in the partition
	numElements = partitionSize;
	
	/* heuristic search looking for better partitions (i.e., partitions that
	   minimize the squared distance to all of the partitions in the
	   list of partitions, prtList) */
	PopPartition *avePart[2];
	avePart[0] = new PopPartition( *partList[0] );
	avePart[1] = new PopPartition( *partList[0] );
	int curAve = 0;
	double bestDist = avePart[0]->ssDistance(partList);
	clock_t startCpuTime = clock();
	float nextTimeGoal = 10.0, goalIncrement = 10.0;
	bool foundBetter = false;
	int passNum = 1;
	bool showedTimeStatus = false;
	cout << "   Searching for mean partition" << endl << endl;
	do
		{
		foundBetter = false;
		for (int i=0; i<numElements; i++)
			{
			int newAve = flip(curAve);
			int deg = avePart[newAve]->getDegree();
			for (int j=1; j<=deg+1; j++)
				{
				avePart[newAve]->setElement(i, j);
				avePart[newAve]->setRgf();
				avePart[newAve]->setDegree();
				avePart[newAve]->setSubsets();
				double d = avePart[newAve]->ssDistance(partList);
				if ( d < bestDist )
					{
					curAve = newAve;
					newAve = flip(curAve);
					bestDist = d;
					foundBetter = true;
					}
				*avePart[newAve] = *avePart[curAve];
				clock_t currentCpuTime = clock();
				float numSecsElapsed = (currentCpuTime - startCpuTime) / (float)CLOCKS_PER_SEC;
				if (numSecsElapsed >= nextTimeGoal)
					{
					string timeStr = MyString::formatTime(numSecsElapsed);
					char tempC[50];
					sprintf(tempC, "%d/%d", i, numElements);
					string progStr = tempC;
					if (showedTimeStatus == false)
						{
						cout << "   " << setw(10) << "Time" << setw(10) << "Pass" << setw(10) << "Progress" << setw(10) << "Score" << endl;
						cout << "   ----------------------------------------" << endl;
						}
					cout << "   " << setw(10) << timeStr << setw(10) << passNum << setw(10) << progStr << fixed << setprecision(2) << setw(10) << bestDist << endl;
					nextTimeGoal += goalIncrement;
					showedTimeStatus = true;
					}
				}
#			if defined (MAC_GUI)
			if ( [appControllerPtr cancelAction] == YES )
				{
				cancelMeanPartitionProcess = true;
				break;
				}
#			endif
			}
		passNum++;
		} while (foundBetter == true);
	if (showedTimeStatus == true)
		cout << "   ----------------------------------------" << endl << endl;
	score = bestDist;

	// allocate the partition
	allocatePartition();
	
	// set up the partition (we set the partition to the same as the first in the list
	for (int i=0; i<numElements; i++)
		part[i] = avePart[curAve]->part[i];

	// set the degree
	degree = avePart[curAve]->degree;

	// put partition into RGF form
	setRgf();

	// get the degree
	setDegree();

	// set the subsets
	setSubsets();
	
	// calculate a vector containing the variance components
	for (int i=0; i<numElements; i++)
		var.push_back( bestDist - ssDistance( partList, i ) );
	double sum = 0.0;
	for (vector<double>::iterator p=var.begin(); p != var.end(); p++)
		sum += (*p);
	for (vector<double>::iterator p=var.begin(); p != var.end(); p++)
		{
		if (sum < 0.0000001)
			(*p) = 0.0;
		else
			(*p) /= sum;
		}
		
	// free memory 
	delete avePart[0];
	delete avePart[1];
	
	// END serial calculation of mean partition
	
#	else

	// BEGIN parallel (threaded) calculation of mean partition

	// check that all of the partitions in the list have the same size
	int partitionSize = partList[0]->getNumElements();
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		if ( (*p)->getNumElements() != partitionSize )
			{
			cerr << "ERROR: PopPartition list has an inconsistent number of elements" << endl;
			exit(1);
			}
		}
		
	// set the number of elements in the partition
	numElements = partitionSize;
	
	// determine the number of processors that are available
	int numProcessors = StMachineInfo::numProcs();
	int numThreads = numProcessors - 2;
	if (numThreads < 2)
		numThreads = 2;
	std::cout << "numProcessors = " << numProcessors << std::endl;
	std::vector<std::vector<int> > elemBlock(numElements/numThreads + 1);
	for (int i=0, j=0; i<numElements; i++)
		{
		elemBlock[j].push_back(i);
		if (elemBlock[j].size() >= numThreads)
			j++;
		}
	/*for (int i=0; i<elemBlock.size(); i++)
		{
		for (int j=0; j<elemBlock[i].size(); j++)
			std::cout << elemBlock[i][j] << " ";
		std::cout << std::endl;
		}*/
	
	// allocate a mean partition
	PopPartition *avePart = new PopPartition( *partList[0] );
	double bestDist = avePart->ssDistance(partList);
	
	// allocate memory for threads
	PopPartition** testParts = new PopPartition*[numElements];
	for (int i=0; i<numElements; i++)
		testParts[i] = new PopPartition( *partList[0] );
	ThreadSettings *threadSettings = new ThreadSettings[numElements];
	int* newElemIndex = new int[numElements];
	double* dists = new double[numElements];

	bool foundBetter = false;
	int passNum = 1;
	bool showedTimeStatus = false, cancelledAnalysis = false;
	cout << "   Searching for mean partition" << endl << endl;
	do
		{
		foundBetter = false;

		for (int i=0; i<elemBlock.size(); i++)
			{
			// spawn a thread for each element to check
			pthread_t* elemThreads = new pthread_t[elemBlock[i].size()];
			for (int j=0; j<elemBlock[i].size(); j++)
				{
				threadSettings[j].elem         = elemBlock[i][j];
				threadSettings[j].partPtr      = testParts[ elemBlock[i][j] ];
				threadSettings[j].list         = &partList;
				threadSettings[j].curBest      = bestDist;
				threadSettings[j].elemIndexPtr = &newElemIndex[ elemBlock[i][j] ];
				threadSettings[j].distPtr      = &dists[ elemBlock[i][j] ];
				pthread_create( &elemThreads[j], NULL, threadDist, (void*) &threadSettings[j] );
				}
			for (int j=0; j<elemBlock[i].size(); j++)
				pthread_join( elemThreads[j], NULL );
			delete [] elemThreads;	
				
#			if defined (MAC_GUI)
			if ( [appControllerPtr cancelAction] == YES )
				{
				cancelledAnalysis = true;
				break;
				}
#			endif
			
			// modify the average partition if better element settings were found
			for (int j=0; j<elemBlock[i].size(); j++)
				{
				if (newElemIndex[ elemBlock[i][j] ] != -1)
					{
					foundBetter = true;
					avePart->setElement( elemBlock[i][j], newElemIndex[ elemBlock[i][j] ] );
					}
				}
			if (foundBetter == true)
				{
				avePart->setRgf();
				avePart->setDegree();
				avePart->setSubsets();
				double d = avePart->ssDistance(partList);
				if (d < bestDist)
					bestDist = d;
				}
#			if defined (MAC_GUI)
			if ( [appControllerPtr cancelAction] == YES )
				{
				cancelledAnalysis = true;
				break;
				}
#			endif
			}
			
		if (cancelledAnalysis == true)
			break;
		for (int i=0; i<numElements; i++)
			std::cout << i << " -- " << newElemIndex[i] << " " << dists[i] << std::endl;
		
		passNum++;
		} while (foundBetter == true);

	if (cancelledAnalysis == false)
		{
		
		bestDist = avePart->ssDistance(partList);
			
		if (showedTimeStatus == true)
			cout << "   ----------------------------------------" << endl << endl;
		score = bestDist;

		// allocate the partition
		allocatePartition();
		
		// set up the partition (we set the partition to the same as the first in the list
		for (int i=0; i<numElements; i++)
			part[i] = avePart->part[i];

		// set the degree
		degree = avePart->degree;

		// put partition into RGF form
		setRgf();

		// get the degree
		setDegree();

		// set the subsets
		setSubsets();
		
		// calculate a vector containing the variance components
		for (int i=0; i<numElements; i++)
			var.push_back( bestDist - ssDistance( partList, i ) );
		double sum = 0.0;
		for (vector<double>::iterator p=var.begin(); p != var.end(); p++)
			sum += (*p);
		for (vector<double>::iterator p=var.begin(); p != var.end(); p++)
			{
			if (sum < 0.0000001)
				(*p) = 0.0;
			else
				(*p) /= sum;
			}
		
		}
			
	// free memory
	delete avePart;
	for (int i=0; i<numElements; i++)
		delete testParts[i];
	delete [] testParts;
	delete [] threadSettings;
	delete [] newElemIndex;
	delete [] dists;

	// END parallel (threaded) calculation of mean partition

#	endif
}

#ifdef THREADED_MEAN_PART

void *threadDist(void *mstPtr) {

	// get settings
	ThreadSettings* mySettings          = (ThreadSettings*)mstPtr;
	int i                               = mySettings->elem;
	PopPartition* part                  = mySettings->partPtr;
	std::vector<PopPartition*> partList = *(mySettings->list);
	double bestDist                     = mySettings->curBest;
	int* elemIndexPtr                   = mySettings->elemIndexPtr;
	double* dist                        = mySettings->distPtr;
	
	std::cout << "Starting thread for element " << i << std::endl;
	
	// make a copy of the partition
	PopPartition* avePart = new PopPartition( *part );
	
	int deg = part->getDegree();
	(*elemIndexPtr) = -1;
	for (int j=1; j<=deg+1; j++)
		{
		*avePart = *part;
		avePart->setElement(i, j);
		avePart->setRgf();
		avePart->setDegree();
		avePart->setSubsets();
		double d = avePart->ssDistance(partList);
		if (j == 1)
			(*dist) = d;
		if ( d < bestDist )
			{
			bestDist = d;
			(*elemIndexPtr) = j;
			(*dist) = d;
			}
#		if defined (MAC_GUI)
		if ( [appControllerPtr cancelAction] == YES )
			break;
#		endif
		}
		
	// free memory
	delete avePart;
	std::cout << "Finishing thread for element " << i << std::endl;
	
	return NULL;
}

#endif

PopPartition::~PopPartition(void) {

	delete [] part;
	deleteSubsets();
}

PopPartition &PopPartition::operator=(PopPartition &p) {

	if (this != &p)
		{
		if (numElements != p.numElements)
			{
			delete [] part;
			part = new int[p.numElements];
			}
		numElements = p.numElements;
		degree = p.degree;
		for (int i=0; i<numElements; i++)
			part[i] = p.part[i];
		setSubsets();
		}
	return *this;
}

bool PopPartition::operator==(const PopPartition &p) {

	if (numElements != p.numElements)
		return false;
	if (degree != p.degree)
		return false;
	for (int i=0; i<numElements; i++)
		if (part[i] != p.part[i])
			return false;
	return true;
}

void PopPartition::allocatePartition(void) {

	part = new int[numElements];
	for (int i=0; i<numElements; i++)
		part[i] = 0;
}

double PopPartition::ssDistance(vector<PopPartition *> &partList) {

	// find the partition with the largest degree
	int largestDegree = getDegree();
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		if ( (*p)->getDegree() > largestDegree )
			largestDegree = (*p)->getDegree();
		}
		
	// allocate a matrix large enough to hold any cost matrix in the list
	int **r = new int*[largestDegree];
	r[0] = new int[largestDegree * largestDegree];
	for (int i=1; i<largestDegree; i++)
		r[i] = r[i-1] + largestDegree;
		
	// allocate a bit field for testing whether the bits in two subsets are both on
	MbBitfield *test = new MbBitfield(numElements);
			
	// calculate the average distance from this partition to all of the paritions in the list
	double distanceSum = 0.0;
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		// get pointers to the two partitions
		PopPartition *part1, *part2;
		if ( this->getDegree() > (*p)->getDegree() )
			{
			part1 = this;
			part2 = (*p);
			}
		else
			{
			part1 = (*p);
			part2 = this;
			}

		// set the number of rows/columns of the cost matrix
		int ncols = part1->getDegree();
		int nrows = part2->getDegree();
		
		// zero-out the cost matrix
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncols; j++)
				r[i][j] = 0;

		// initialize the cost matrix
		for (int i=0; i<nrows; i++)
			{
			MbBitfield *bf2 = part2->getSubset(i);
			for (int j=0; j<ncols; j++)
				{
				MbBitfield *bf1 = part1->getSubset(j);
				(*test) = (*bf1) & (*bf2);
				for (int k=0; k<numElements; k++)
					{
					if ( test->isBitSet(k) == true )
						r[i][j]++;
					}
				}
			}

		// calculate assignment cost
		Hungarian *h = new Hungarian(r, nrows, ncols);
		int d = numElements - h->getAssignmentCost();
		
		// add to the sum
		distanceSum += (d * d);
		
		delete h;
		}
		
	// free memory
	delete [] r[0];
	delete [] r;
	delete test;
	
	// return the average distance
	return distanceSum / partList.size();
}

double PopPartition::ssDistance(vector<PopPartition *> &partList, int elemToDelete) {

	// make a list of partitions with one of the elements deleted
	vector<PopPartition *> tempPartList;
	for (vector<PopPartition *>::iterator p=partList.begin(); p != partList.end(); p++)
		{
		PopPartition *newPart = new PopPartition(*(*p));
		newPart->deleteElement(elemToDelete);
		tempPartList.push_back(newPart);
		//cout << "Pair " << elemToDelete << " " << newPart->getNumElements() << endl;
		//(*p)->print();
		//newPart->print();
		}
		
	// make a copy of myself
	PopPartition *selfPart = new PopPartition(*this);
	selfPart->deleteElement(elemToDelete);
		
	// find the partition with the largest degree
	int largestDegree = selfPart->getDegree();
	for (vector<PopPartition *>::iterator p=tempPartList.begin(); p != tempPartList.end(); p++)
		{
		if ( (*p)->getDegree() > largestDegree )
			largestDegree = (*p)->getDegree();
		}
		
	// allocate a matrix large enough to hold any cost matrix in the list
	int **r = new int*[largestDegree];
	r[0] = new int[largestDegree * largestDegree];
	for (int i=1; i<largestDegree; i++)
		r[i] = r[i-1] + largestDegree;
		
	// allocate a bit field for testing whether the bits in two subsets are both on
	MbBitfield *test = new MbBitfield(selfPart->getNumElements());
			
	// calculate the average distance from this partition to all of the paritions in the list
	double distanceSum = 0.0;
	for (vector<PopPartition *>::iterator p=tempPartList.begin(); p != tempPartList.end(); p++)
		{
		// get pointers to the two partitions
		PopPartition *part1, *part2;
		if ( selfPart->getDegree() > (*p)->getDegree() )
			{
			part1 = selfPart;
			part2 = (*p);
			}
		else
			{
			part1 = (*p);
			part2 = selfPart;
			}

		// set the number of rows/columns of the cost matrix
		int ncols = part1->getDegree();
		int nrows = part2->getDegree();
		
		// zero-out the cost matrix
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncols; j++)
				r[i][j] = 0;

		// initialize the cost matrix
		for (int i=0; i<nrows; i++)
			{
			MbBitfield *bf2 = part2->getSubset(i);
			for (int j=0; j<ncols; j++)
				{
				MbBitfield *bf1 = part1->getSubset(j);
				(*test) = (*bf1) & (*bf2);
				for (int k=0; k<part1->getNumElements(); k++)
					{
					if ( test->isBitSet(k) == true )
						r[i][j]++;
					}
				}
			}

		// calculate assignment cost
		Hungarian *h = new Hungarian(r, nrows, ncols);
		int d = part1->getNumElements() - h->getAssignmentCost();
		
		// add to the sum
		distanceSum += (d * d);
		
		delete h;
		}

	// calculate the average distance
	double aveDistance = distanceSum / tempPartList.size();
		
	// free memory
	delete [] r[0];
	delete [] r;
	delete test;
	for (vector<PopPartition *>::iterator p=tempPartList.begin(); p != tempPartList.end(); p++)
		delete (*p);
		
	// return the average distance
	return aveDistance;
}

void PopPartition::deleteElement(int elemToDelete) {

	int *tempPart = new int[numElements - 1];
	for (int i=0, j=0; i<numElements; i++)
		{
		if (i != elemToDelete)
			tempPart[j++] = part[i];
		}
	delete [] part;
	part = tempPart;
	numElements--;
	setRgf();
	setDegree();
	deleteSubsets();
	setSubsets();
}

void PopPartition::deleteSubsets(void) {

	for (vector<MbBitfield *>::iterator b=subsets.begin(); b != subsets.end(); b++)
		delete (*b);
	subsets.clear();
}

int PopPartition::distance(PopPartition *p) {

	PopPartition *p1 = this;
	PopPartition *p2 = p;
	
	if (p1 == p2)
		return 0;
	
	if ( p1->getNumElements() != p2->getNumElements() )
		return -1;
	
	// partition 1 should have a larger degree than partition 2
	PopPartition *part1, *part2;
	if ( p1->getDegree() > p2->getDegree() )
		{
		part1 = p1;
		part2 = p2;
		}
	else
		{
		part1 = p2;
		part2 = p1;
		}
		
	// set the number of rows/columns of the cost matrix
	int ncols = part1->getDegree();
	int nrows = part2->getDegree();
		
	// allcoate the cost matrix
	int **r = new int*[nrows];
	r[0] = new int[nrows * ncols];
	for (int i=1; i<nrows; i++)
		r[i] = r[i-1] + ncols;
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			r[i][j] = 0;
			
	// initialize the cost matrix
	MbBitfield *test = new MbBitfield(numElements);
	for (int i=0; i<nrows; i++)
		{
		MbBitfield *bf2 = part2->getSubset(i);
		for (int j=0; j<ncols; j++)
			{
			MbBitfield *bf1 = part1->getSubset(j);
			(*test) = (*bf1) & (*bf2);
			for (int k=0; k<numElements; k++)
				{
				if ( test->isBitSet(k) == true )
					r[i][j]++;
				}
			}
		}
	delete test;
		
#	if 0
	for (int i=0; i<nrows; i++)
		{
		for (int j=0; j<ncols; j++)
			std::cout << r[i][j] << ",";
		std::cout << std::endl;
		}
	std::cout << "done" << std::endl;
	getchar();
#	endif

	// calculate assignment cost
	Hungarian *h = new Hungarian(r, nrows, ncols);
	int d = numElements - h->getAssignmentCost();
	
	// free memory
	delete [] r[0];
	delete [] r;
	delete h;

	return d;
}

int PopPartition::flip(int x) {

	if (x == 0)
		return 1;
	return 0;
}

void PopPartition::print(void) {

	
	for (int i=0; i<numElements; i++)
		{
		//stprint("%d", part[i] );
		//if (i < numElements - 1)
		//	stprint(",");
		cout << part[i];
		if (i < numElements - 1)
			cout << ",";
		}
	//stprint("\n");
	cout << endl;
}

void PopPartition::setSubsets(void) {

	if ( subsets.size() > 0 )
		deleteSubsets();

	for (int i=0; i<degree; i++)
		{
		MbBitfield *bf = new MbBitfield(numElements);
		int idx = i+1;
		for (int j=0; j<numElements; j++)
			{
			if ( part[j] == idx )
				bf->setBit(j);
			}
		subsets.push_back( bf );
		}
}

void PopPartition::setDegree(void) {

	/* Find the largest number in the partition. This number will be the
	   degree of the model. This calculation assumes that the partition
	   is in the restricted growth function (RGF) format. */

	int largestNumFound = 0;
	int smallestNumFound = part[0];
	for (int i=0; i<numElements; i++)
		{
		if (part[i] > largestNumFound)
			largestNumFound = part[i];
		}
	degree = largestNumFound - smallestNumFound + 1;
}

void PopPartition::setRgf(void) {

	int *rgf = new int[numElements];
	for (int i=0; i<numElements; i++)
		rgf[i] = -1;
	
	int idx = 1;
	for (int i=0; i<numElements; i++)
		{
		if (rgf[i] == -1)
			{
			int partId = part[i];
			for (int j=0; j<numElements; j++)
				{
				if (part[j] == partId)
					rgf[j] = idx;
				}
			idx++;
			}
		}
	
	for (int i=0; i<numElements; i++)
		{
		if (rgf[i] == -1)
			{
			cerr << "ERROR: Problem initializing partition" << endl;
			exit(1);
			}
		part[i] = rgf[i];
		}
	
	delete [] rgf;
}




