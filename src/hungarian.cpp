#include "hungarian.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

#define HUNGARIAN_NOT_ASSIGNED			0 
#define HUNGARIAN_ASSIGNED 				1
#define INF								(0x7FFFFFFF)

using namespace std;



Hungarian::Hungarian(int **m, int nr, int nc) {

	/* allocate minimizedMatrix */
	minimizedMatrix = new int*[ imax(nr, nc) ];
	for (int i=0; i<imax(nr, nc); i++)
		minimizedMatrix[i] = new int[ imax(nr, nc) ];

	/* initialize the hungarian problem using the cost matrix */
	matrixSize = hungarianInit(m, nr, nc);

	/* solve the assignement problem */
	hungarianSolve();

	/* calculate the cost */
	assignmentCost = 0;
	for (int i=0; i<matrixSize; i++)
		{
		for (int j=0; j<matrixSize; j++)
			if ( assignment[i][j] == 1 )
				assignmentCost += minimizedMatrix[i][j];
		}

	/* free memory */
	delete [] cost[0];
	delete [] cost;
	delete [] assignment[0];
	delete [] assignment;
	for (int i=0; i<matrixSize; i++)
		delete [] minimizedMatrix[i];
	delete [] minimizedMatrix;
}

int Hungarian::hungarianInit(int **costMatrix, int rows, int cols) {

	int maxCost = 0;
	int orgCols = cols;
	int orgRows = rows;

	/* is the number of cols not equal to number of rows : if yes, expand with 0-cols / 0-cols */
	rows = imax(cols, rows);
	cols = rows;

	numRows = rows;
	numCols = cols;

	cost = new int*[rows];
	cost[0] = new int[rows * cols];
	for (int i=1; i<rows; i++)
		cost[i] = cost[i-1] + cols;
	assignment = new int*[rows];
	assignment[0] = new int[rows * cols];
	for (int i=1; i<rows; i++)
		assignment[i] = assignment[i-1] + cols;

	for (int i=0; i<numRows; i++) 
		{
		for (int j=0; j<numCols; j++) 
			{
			cost[i][j] =  (i < orgRows && j < orgCols) ? costMatrix[i][j] : 0;
			assignment[i][j] = 0;
			if (maxCost < cost[i][j])
				maxCost = cost[i][j];
			minimizedMatrix[i][j] = cost[i][j];
			}
		}

	for(int i=0; i<numRows; i++) 
		for(int j=0; j<numCols; j++) 
			cost[i][j] =  maxCost - cost[i][j];

	return rows;
}

void Hungarian::hungarianSolve(void) {

	int k, l, s, t, q;
	
	int myCost = 0;
	int m = numRows;
	int n = numCols;

	int *colMate     = new int[4*numRows + 4*numCols];
	int *unchosenRow = colMate     + numRows;
	int *rowDec      = unchosenRow + numRows;
	int *slackRow    = rowDec      + numRows;
	int *rowMate     = slackRow    + numRows;
	int *parentRow   = rowMate     + numCols;
	int *colInc      = parentRow   + numCols;
	int *slack       = colInc      + numCols;

	for (int i=0; i<numRows; i++) 
		{
		colMate[i]     = 0;
		unchosenRow[i] = 0;
		rowDec[i]      = 0;
		slackRow[i]    = 0;
		}
	for (int j=0; j<numCols; j++) 
		{
		rowMate[j]   = 0;
		parentRow[j] = 0;
		colInc[j]    = 0;
		slack[j]     = 0;
		}

	for (int i=0; i<numRows; ++i)
		for (int j=0; j<numCols; ++j)
			assignment[i][j] = HUNGARIAN_NOT_ASSIGNED;

	for (l=0; l<n; l++)
		{
		s = cost[0][l];
		for (k=1; k<m; k++) 
			if ( cost[k][l] < s )
				s = cost[k][l];
		myCost += s;
		if ( s != 0 )
			for (k=0; k<m; k++)
				cost[k][l] -= s;
		}

	t = 0;
	for (l=0; l<n; l++)
		{
		rowMate[l] = -1;
		parentRow[l] = -1;
		colInc[l] = 0;
		slack[l] = INF;
		}
	for (k=0; k<m; k++)
		{
		s = cost[k][0];
		for (l=1; l<n; l++)
			if (cost[k][l] < s)
				s = cost[k][l];
		rowDec[k] = s;
		for (l=0; l<n; l++)
			{
			if ( s == cost[k][l] && rowMate[l] < 0 )
				{
				colMate[k] = l;
				rowMate[l] = k;
				goto rowDone;
				}
			}
		colMate[k] = -1;
		unchosenRow[t++] = k;
		rowDone:
		;
		}

	int unmatched = t;
	if (t == 0)
		goto done;
	while (1)
		{
		q = 0;
		while (1)
			{
			while ( q < t )
				{
				{
				k = unchosenRow[q];
				s = rowDec[k];
				for (l=0; l<n; l++)
					{
					if (slack[l])
						{
						int del;
						del = cost[k][l] - s + colInc[l];
						if ( del<slack[l] )
							{
							if (del == 0)
								{
								if (rowMate[l] < 0)
									goto breakthru;
								slack[l] = 0;
								parentRow[l] = k;
								unchosenRow[t++] = rowMate[l];
								}
							else
								{
								slack[l] = del;
								slackRow[l] = k;
								}
							}
						}
					}
				}
			q++;
			}

			s = INF;
			for (l=0; l<n; l++)
				if (slack[l] && slack[l] < s)
					s = slack[l];
			for (q=0; q<t; q++)
				rowDec[ unchosenRow[q] ] += s;
			for (l=0; l<n; l++)
				{
				if (slack[l])
					{
					slack[l] -= s;
					if (slack[l] == 0)
						{
						k = slackRow[l];
						if (rowMate[l] < 0)
							{
							for (int j=l+1; j<n; j++)
								if (slack[j] == 0)
									colInc[j] += s;
							goto breakthru;
							}
						else
							{
							parentRow[l] = k;
							unchosenRow[t++] = rowMate[l];
							}
						}
					}
				else
					colInc[l] += s;
				}
			}
		breakthru:
		while (1)
			{
			int j = colMate[k];
			colMate[k] = l;
			rowMate[l] = k;
			if ( j < 0 )
				break;
			k = parentRow[j];
			l = j;
			}
		if ( --unmatched == 0 )
			goto done;
		t = 0;
		for (l=0; l<n; l++)
			{
			parentRow[l] = -1;
			slack[l] = INF;
			}
		for (k=0; k<m; k++)
			if (colMate[k] < 0)
				unchosenRow[t++] = k;
		}
	done:

	for (k=0; k<m; k++)
		for (l=0; l<n; l++)
			if ( cost[k][l] < rowDec[k] - colInc[l] ) {
				cout << "help 1 " << endl;
				getchar();
				exit(0);
			}
	for (k=0; k<m; k++) {
		l = colMate[k];
		if ( l < 0 || cost[k][l] != rowDec[k] - colInc[l] ) {
			cout << "help 2 " << endl;
			getchar();
			exit(0);
	  }
	}
	k = 0;
	for (l=0; l<n; l++)
		if (colInc[l])
			k++;
	if (k > m)
		{
		cout << "help 3 " << endl;
		getchar();
		exit(0);
		}

	for (int i=0; i<m; ++i)
		assignment[i][ colMate[i] ] = HUNGARIAN_ASSIGNED;
	for (k=0; k<m; ++k)
		for (l=0; l<n; ++l)
			cost[k][l] = cost[k][l] - rowDec[k] + colInc[l];
	for (int i=0; i<m; i++)
		myCost += rowDec[i];
	for (int i=0; i<n; i++)
		myCost -= colInc[i];
	
	delete [] colMate;
}




