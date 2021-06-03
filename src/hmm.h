#ifndef HMM_H_INCLUDED
#define HMM_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#define N 15
#define M 2
#define L 10
#define C 1
#define T 800000
#define RANDMAX 1000


/**Function declaration**/
int printobservation(int arr[][C], int len,int width);
int printvector(double arr[], int len);
int printmatrix(double **mat, int nrow, int ncol);
int initializer(double prob[], int n, int choice);
int genobslen(int numt[], int numc, int totallen, int choice);
int balance(int arr[], int first, int last, int total);
double scaler(double *arr, int len);
double forwardstep(double** alpha, double a[][N], double b[][L][M], double pi[], int** o, int n, int numt, int choice);
double backwardstep(double** beta, double a[][N], double b[][L][M], double pi[], int** o, int n, int numt, int choice);
int getgamma(double** gamma, double** alpha, double** beta, int n, int numt);
int gethiddenstates(int hiddeni[],double **gamma,int numt,int n);
int getxi(double*** xi, double** alpha, double** beta, double a[][N], double b[][L][M], int** o, int n, int numt);
int getnumeratordenominator(double pi_numerator[], double *pi_denominator, double a_numerator[][N],
                            double a_denominator[], double b_numerator [][L][M], double b_denominator[],
                            int** o, double** gamma, double*** xi,int n,int m,int numt);
int hmm(int n, int m, int numc, int numt[C], int*** o,
		double a[N][N], double b[N][L][M], double pi[N],
		double*** alpha, double*** beta, double*** gamma,
		double**** xi, double pi_numerator[C][N], double pi_denominator[C],
		double a_numerator[C][N][N], double a_denominator[C][N],
		double b_numerator[C][N][L][M], double b_denominator[C][N],
		double numerator, double denominator,
		int sigmat, int par, int tn, unsigned int seed,
		int choicegetobs, int choiceinitial, int choicefb, clock_t timearr[],
		int *piter, double *plo, double *pl,
		int maxiter, double eps, int endwithouteps);

#endif // HMM_H_INCLUDED
