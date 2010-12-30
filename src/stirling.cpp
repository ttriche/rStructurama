#include "integer.h"
#include "stirling.h"
#include "st.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace std;



Stirling::Stirling(int n) {

	num = n;
	
	/* allocate memory */
	Integer **s1[2];
	s1[0] = new Integer*[2*(n+1)];
	s1[1] = s1[0] + (n+1);
	for (int i=0; i<2; i++)
		for (int j=0; j<n+1; j++)
			s1[i][j] = NULL;
	int curCol = 0, prevCol = 1;
				
	lnStirling = new double[ num + 1 ];
	for (int i=0; i<num+1; i++)
		lnStirling[i] = 0.0;

	priorOnK = new double [num];
	for (int i=0; i<num; i++)
		priorOnK[i] = 0.0;
		
	/* calculate stirling numbers of the first kind */
	clock_t startCpuTime = clock();
	float nextTimeGoal = 20.0, goalIncrement = 10.0;
	bool showedTimeStatus = false;
	for (int i=1, k=0; i<=n; i++)
		{
		for (int j=1; j<=i; j++)
			{
			if (s1[curCol][j] != NULL)
				delete s1[curCol][j];
			s1[curCol][j] = new Integer();
			
			//cout << i << " " << j << endl;
			if (j == 1)
				{
				/* case 1 */
				if (i == 1)
					{
					/*s1[i][j] = 1;*/
					(*s1[curCol][j]) = 1;
					}
				else
					{
					/*s1[i][j] = s1[i-1][j] * (i-1);*/
					(*s1[curCol][j]) = (*s1[prevCol][j]) * (i-1);
					}
				}
			else if (i == j)
				{
				/* case 2 */
				/*s1[i][j] = 1;*/
				(*s1[curCol][j]) = 1;
				}
			else if (j+1 == i)	
				{
				/* case 3 */
				/*s1[i][j] = i * (i-1) / 2;*/
				(*s1[curCol][j]) = (i * (i-1) / 2);
				}
			else
				{
				/* case 4 */
				/*s1[i][j] = s1[i-1][j-1] + (i-1) * s1[i-1][j];*/
				Integer *y = new Integer();
				(*y) = (*s1[prevCol][j]) * (i-1);
				(*s1[curCol][j]) = (*s1[prevCol][j-1]) + (*y);
				delete y;
				}

			/* calculate the log stirling number */
			if (i == n)
				{
				lnStirling[j] = (*s1[curCol][j]).ln();
				//cout << "S1(" << j << ") = " << (*s1[curCol][j]) << " " << lnStirling[j] << endl;
				}
				
			// print status
			clock_t currentCpuTime = clock();
			float numSecsElapsed = (currentCpuTime - startCpuTime) / (float)CLOCKS_PER_SEC;
			k++;
			if (numSecsElapsed >= nextTimeGoal)
				{
				string timeStr = MyString::formatTime(numSecsElapsed);
				char tempC[50];
				sprintf(tempC, "%d/%d", k, (n+1)*n/2);
				string progStr = tempC;
				if (showedTimeStatus == false)
					{
					cout << "   " << setw(10) << "Time" << setw(10) << "Progress" << endl;
					cout << "   --------------------" << endl;
					}
				cout << "   " << setw(10) << timeStr << setw(10) << progStr << endl;
				nextTimeGoal += goalIncrement;
				showedTimeStatus = true;
				}
			
			}
			
		if (curCol == 0)
			{
			curCol = 1;
			prevCol = 0;
			}
		else
			{
			curCol = 0;
			prevCol = 1;
			}
			
		}
	if (showedTimeStatus == true)
		cout << "   --------------------" << endl;

	/* free memory */
	for (int i=0; i<2; i++)
		for (int j=0; j<n+1; j++)
			if (s1[i][j] != NULL)
				delete s1[i][j];
	delete [] s1[0];
}

Stirling::~Stirling(void) {

	delete [] lnStirling;
	lnStirling = NULL;
	delete [] priorOnK;
}

void Stirling::calculateProbK(double alpha) {

	double ax = 0.0;
	for (int i=1; i<=num; i++)
		ax += log(alpha + i - 1.0);
		
	for (int i=0; i<num; i++)
		priorOnK[i] = exp(lnStirling[i+1] + (i+1) * log(alpha) - ax);
}

void Stirling::calculateLnProbK(double alpha) {

	double ax = 0.0;
	for (int i=1; i<=num; i++)
		ax += log(alpha + i - 1.0);
		
	for (int i=0; i<num; i++)
		priorOnK[i] = lnStirling[i+1] + (i+1) * log(alpha) - ax;
}

