#include "hmm.h"

//int getxi(double xi[][N][N], double alpha[][N], double beta[][N], double a[][N], double b[][L][M], int** o, int n, int numt)
int getxi(double*** xi, double** alpha, double** beta, double a[][N], double b[][L][M], int** o, int n, int numt)
{
    int i,j,t;
    double sum;
   //double x[N][T];
	double** x = (double**)malloc(N * sizeof(double*));
	for (int i0 = 0; i0 < N; i0++) {
		x[i0] = (double*)calloc(T, sizeof(double));
	}
    for(i=0; i<n; i++)
        for(t=0; t<numt; t++)
            for(int l=0; l<L; l++)
            {
                x[i][t]=1;
                x[i][t]*=b[i][l][o[t][l]];
            }
    for(t=0; t<numt-1; t++)
    {
        for(i=0,sum=0; i<n; i++)
            for(j=0; j<n; j++)
            {
                xi[t][i][j]=alpha[t][i]*a[i][j]*x[j][t+1]*beta[t+1][j];
                sum+=xi[t][i][j];
            }
        for(i=0; i<n; i++)
            for(j=0; j<n; j++)
                xi[t][i][j]/=sum;
    }
	for (int i0 = 0; i0 < N; i0++)
		free(x[i0]);
	free(x);
    return(0);
}
