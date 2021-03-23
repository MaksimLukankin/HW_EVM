#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include "headers.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <sys/sysinfo.h>
using namespace std;
void vvod_b(double* b, double* mat, int n, int first_row, int last_row, int p)
   {
        double a = 0;
        for(int i = first_row; i <= last_row; i++)
		{
            for(int j = 0; j < (n-1) / 2 + 1; j++)
			{
                a = a + mat[n * (i - first_row) + 2 * j];
            }
            b[i - first_row] = a;
            a = 0;
        }	
   }

int main(int argc, char * argv[]) 
{
   int p;  // номер потока
   int num; // число потоков
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&num);
   MPI_Comm_rank(MPI_COMM_WORLD,&p);
   int err = 0; 
   int res1 = 0, res2 = 0;
   string filename = ""; 
   int n=0, m=0, k=0;
   int i;
   double start = 0.0, end = 0.0;
   struct sysinfo info;
   int first_row = 0, last_row = 0;
   if (argc == 5) 
   {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
        if(k<0||k>4||n<1||m<1)
        {
           printf("Error. Incorrect arguments");
           MPI_Finalize();
           return 0;
        }
        filename = argv[4];
    }
    else if (argc == 4) 
    {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
       if(k<0||k>4||n<1||m<1)
        {
           cout<<"Error. Incorrect arguments";
           MPI_Finalize();
           return 0;
        }
        if(k == 0)
        {
           printf("Error. No file name");
           MPI_Finalize();
           return 0;
        }     
    }
    else 
    {
        printf("Error. Incorrect arguments");
        MPI_Finalize();
        return 0;
    }
    if (num == 0)
    {
      printf("Error. 0 prosesses");
      MPI_Finalize();
      return 0;
    }
    first_row = n * p / num;
    last_row = n * (p + 1) / num - 1;
    sysinfo (&info);
    long unsigned int ramflag;
    ramflag = info.freeram;
    if(ramflag <(n*sizeof(int) + n * (last_row - first_row + 5)*sizeof(double))) 
    {
       cout << ramflag << "\n";
       cout << "Error. Not enough ram";
       MPI_Finalize();
       return 0;
    }
    double* mat = new double[n * (last_row - first_row + 1)];
    double* b = new double[last_row - first_row + 1];
    double* x = new double[n];
    double* per1 = new double[n + 1];
    double* per2 = new double[n + 1];
    per1[0] = 0;
    per2[0] = 0;
    int* ind = new int[n];
    if (k == 0) 
	{
	     err = vvod_file(n, mat, filename, first_row, last_row);
	     if (err && p == 0)
		{
		    cout << "Error: bad file "<< endl;
                    free(mat);
	            free(x); 
	            free(b);
	            free(ind);
                    free(per1);
                    free(per2);
		    MPI_Finalize();
                    return 0;
		}		   
	}
     else
	{
	     vvod_formula(n, k, mat, first_row, last_row);
	}
     vvod_b(b, mat, n,first_row, last_row,p);
     for (i = 0; i < n; i++) 
	{
	    ind[i] = i;
            x[i] = 0.0;
	} 
     MPI_Barrier(MPI_COMM_WORLD);
    for (int h = 0; h < num; h++) {
			if (h == p) {
				for (int i = first_row; i <= last_row && i < m; i++) 
                                {
					for (int j = 0; j < m; j++) 
                                        {
					         cout << mat[(i - first_row) * n + j] << " ";
				        }
				         cout << endl;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
     }
    MPI_Barrier(MPI_COMM_WORLD);
    if(p == 0)
                        {
                            cout << endl;
                        }   
    for (int h = 0; h < num; h++) {
			if (h == p) {
				for (int i = first_row; i <= last_row && i < m; i++) 
                                {
					         cout << b[(i - first_row)] << " ";
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
     }
     MPI_Barrier(MPI_COMM_WORLD);
     if(p == 0)
                        {
                            cout << endl;
                        }
     start = MPI_Wtime();
     err = Jordan_mpi(n, mat, b, ind, p, num, per1, per2, &err);
     MPI_Barrier(MPI_COMM_WORLD);
     end = MPI_Wtime();
     if (p == 0) 
	{
	  cout << "Time: " << end - start<<"\n" << endl;
	}
     if (err == -3) 
         {
           if (p == 0) printf("Error. Det is 0");
           free(mat);
	   free(x); 
	   free(b);
	   free(ind);
           free(per1);
           free(per2);
           MPI_Finalize();
           return 0;
         }
      if (k == 0) 
	   {
		 vvod_file(n, mat, filename, first_row, last_row);
	   }
      else
	   {
		vvod_formula(n, k, mat, first_row, last_row);
	   }  
      MPI_Barrier(MPI_COMM_WORLD);
     double a;
      for (int h = 0; h < num; h++)  
     {
	for (int i = 0; i < n; i++) 
              {
                  a = x[i]; 
		  MPI_Bcast(&a, 1, MPI_DOUBLE, h, MPI_COMM_WORLD);
                  x[i] = a;
	      }
        for (int i = first_row; i <= last_row; i++)
              {
                 x[ind[i]] = b[i - first_row];               
              }
			MPI_Barrier(MPI_COMM_WORLD);
     }
     MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     vvod_b(b, mat, n,first_row, last_row,p);
      for (int h = 0; h < num; h++) {
			if (h == p) {
				for (int i = first_row; i <= last_row && i < m; i++) 
                                {
					for (int j = 0; j < m; j++) 
                                        {
					         cout << mat[(i - first_row) * n + j] << " ";
				        }
				         cout << endl;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
     }
    MPI_Barrier(MPI_COMM_WORLD);
    if(p == 0)
                        {
                            cout << endl;
                        }  
   for (int h = 0; h < num; h++) {
			if (h == p) {
				for (int i = first_row; i <= last_row && i < m; i++) 
                                {
					         cout << b[(i - first_row)] << " ";
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
     }
    
     MPI_Barrier(MPI_COMM_WORLD);
     if(p == 0)
                        {
                            cout << endl;
                        }
       if (p == 0) 
	    {
		printf("Answer:\n");
		for(i = 0; i < m; i++)
		{
                   //x[i] = 0;
		   cout << x[i] << " ";
	        }   
		cout << "\n"<< endl;
	    }
       //MPI_Barrier(MPI_COMM_WORLD);
       double res = normaM(mat, b, x, n, p, num,res1, res2);
       if(p == 0)
           {
             cout << "norma neviazki: " << res << "\n"<< endl;
           }
        //MPI_Barrier(MPI_COMM_WORLD);
        if(p == 0)
           {
             cout << "answer norma: " << norma(x, n)<< "\n"<< endl;
           }
        //MPI_Barrier(MPI_COMM_WORLD);
        free(mat);
        free(x); 
        free(b);
        free(ind);
        free(per1);
        free(per2);
        MPI_Finalize();
        return 0;
}
