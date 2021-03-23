#include <fstream>
#include <iostream>
#include "headers.h"
#include <mpi.h>
using namespace std;
double f(int k, int n, int i, int j) 
{
	double a = 0;
	switch(k) 
	{
		case 1:
			a = n - max(i + 1,j + 1) + 1;
			break;
		case 2:
			a = max(i + 1,j + 1);
			break;
		case 3:
			a = abs(i - j);
			break;
		case 4:
			a = (double) 1 / (i + j + 1);
		break;
	}
	return a;
}
int vvod_formula(int n, int k, double* mat, int first_row, int last_row) 
{
	for (int i = first_row; i <= last_row; i++)  
	{
		for (int j = 0; j < n; j++) 
		{
			mat[(i - first_row) * n + j] = f(k, n, i, j);
		}
	}
	return 0;
}

int vvod_file(int n, double *mat, string filename, int first_row, int last_row) 
{
	double d = 0;
	ifstream file(filename);
	if (!file)
	{
		cout << "Error. File cant be open" << endl;
		return -1000;
	}	
	int i = 0;
	while(!file.eof()) 
	{
		file >> d;
		if (first_row * n <= i && i <= last_row * n + n - 1) 
		{
			mat[(i / n - first_row) * n + i % n] = d;
		}
		i++;
	}
    if (n * n != i - 1)
		return n*n - (i - 1);
    return 0;
}

