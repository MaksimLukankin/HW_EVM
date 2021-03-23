#include <iostream>
#include <cmath>
#include <mpi.h>
#include "headers.h"
using namespace std;
double normaM(double* mat, double* b, double* x, int n, int p, int num, double res1,double res2)
{
    int i;
    int first_row = n * p / num, last_row = (n * (p + 1) / num - 1);
    double s = 0.0, y = 0.0, t = 0.0;
    for (i = first_row; i <= last_row; i++)
    {
        for(int j = 0; j < n; j++)
         { 
           if(j == 0) y = 0.0;
           y = y + mat[(i - first_row) * n + j] * x[j];
           if(j == n - 1)
           {
              s = s + (y - b[i - first_row])*(y - b[i - first_row]);
              t = t + b[i - first_row] * b[i - first_row];
           }
         }    
    }
    MPI_Allreduce(&s,&res1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&t,&res2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    res1 = sqrt(res1);
    res2 = sqrt(res2);
    return res1 / res2;
}
double norma(double* x, int n)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++)
    {
        if (i % 2 == 0)
        {
            s = s + (x[i] - 1) * (x[i] - 1);
        }
        if (i % 2 == 1)
        {
            s = s + (x[i]) * (x[i]);
        }
    }
    s = sqrt(s);
    return s;
}
double normaMat(double* mat, int n,int p, int num, double res1)
{
   double max = 0.0, s = 0.0;
   int i;
   int first_row = n * p / num, last_row = (n * (p + 1) / num - 1);
   for (i = first_row; i <= last_row; i++)
    {
        for(int j = 0; j < n; j++)
         {
           s = s + abs(mat[i - first_row]);
           if(j == n - 1)
           {
              if(max < s) max = s;
              s=0;
           }
         }    
    }
   MPI_Allreduce(&max,&res1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   return res1;
}
int Jordan_mpi(int n, double* mat, double* b, int * ind, int p, int num, double * per1 ,  double * per2,int  *err) 
{
	int first_row = n * p / num, last_row = (n * (p + 1) / num - 1);
	int temp[3];
        double res1 = 1;
        double eps = min(1e-12,1e-12 * normaMat(mat,n,p,num,res1));
        double tp;
        int k, j;
        int sign = 0;
	int i = 0;
	while (i < n)
	{
		int c_max = i, r_max = max(i, first_row);
		for (int k = max(i, first_row) ; k <= last_row; k++) 
		{ 
			for (int l = i; l < n; l++) 
			{ 
				if (fabs(mat[(k - first_row) * n + l]) > fabs(mat[(r_max - first_row) * n + c_max])) 
				{
					r_max = k;
					c_max = l;
				}
			}
		}
		struct double_int 
		{
		    double val;
		    int rank;
		} local_max, global_max;
		MPI_Barrier(MPI_COMM_WORLD);
		if (r_max > last_row)
			local_max.val = 0;
		else
			local_max.val = fabs(mat[(r_max - first_row) * n + c_max]);
		local_max.rank = p;
		MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (p == global_max.rank) 
		{
			temp[0] = r_max;
			temp[1] = c_max;
                        if(mat[(r_max - first_row) * n + c_max] > eps)
                           temp[2] = 1;
                        else 
                           temp[2] = -1;
		}
		MPI_Bcast(temp, 3, MPI_INT, global_max.rank, MPI_COMM_WORLD);
		r_max = temp[0];
		c_max = temp[1];
                sign  = temp[2];
		if (global_max.val < eps) 
		{
		    *err = -3;
		    return -3; 
		}
		for (int l = first_row; l <= last_row; l++)
		{
			swap(mat[(l - first_row) * n + i], mat[(l - first_row) * n + c_max]); 
		}
		swap(ind[i], ind[c_max]);
		int rowI = 0; 
		while (n * (rowI + 1) / num - 1 < i) rowI++;
		int rowMax = 0; 
		while (n * (rowMax + 1) / num - 1 < r_max) rowMax++; 
		if (p == rowI) 
		{
			for (int l = 0; l < n; l++) 
			{
				per1[l] = mat[(i - first_row) * n + l];		
			}
			per1[n] = b[i - first_row];
		}
		MPI_Bcast(per1, n + 1, MPI_DOUBLE, rowI, MPI_COMM_WORLD);
		if (p == rowMax) 
		{

			for (int l = 0; l < n; l++) 
			{
				per2[l] = mat[(r_max - first_row) * n + l];		
			}
			per2[n] = b[r_max - first_row];
			for (int l = 0; l < n; l++) 
			{
				mat[(r_max - first_row) * n + l] = per1[l];	
			}
			b[r_max - first_row] = per1[n];
		}
		MPI_Bcast(per2,n + 1, MPI_DOUBLE, rowMax, MPI_COMM_WORLD);
		if (p == rowI) 
		{
			for (int l = 0; l < n; l++) 
			{
				mat[(i - first_row) * n + l] = per2[l] / global_max.val * sign;
			}
			b[i - first_row] = per2[n] / global_max.val * sign;
		}
		for (j = first_row; j <= last_row; j++) 
                {
		    if(j != i && fabs(mat[(j - first_row) * n + i]) > eps)
                   { 
                    tp = mat[(j - first_row) * n + i];
		    for (k = i; k < n; k++) 
		     {
			mat[(j - first_row) * n + k] -= (per2[k] / global_max.val * sign)  * tp;
		     }
                     b[j - first_row] = b[j - first_row] - per2[n] / global_max.val * sign * tp;	 
                    }               
		}
		i++;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}

