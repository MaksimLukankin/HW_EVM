#include <string>
using namespace std;
int Jordan_mpi(int n, double* mat, double* b, int * ind, int p, int num, double * per1 ,  double * per2,int  *err);
//void true_ind(double* x, double* b, int* ind, int n, int first_row, int last_row, int p);
double f(int k, int n, int i, int j);
int vvod_formula(int n, int k, double* mat, int first_row, int last_row);
int vvod_file(int n, double *mat, string filename, int first_row, int last_row);
void vvod_b(double* b, double* mat, int n, int first_row, int last_row);
double normaM(double* mat, double* b, double* x, int n, int p, int num, double res1,double res2);
double norma(double* x, int n);
double normaMat(double* mat, int n,int p, int num, double res1);



