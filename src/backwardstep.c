#include "hmm.h"

/**Return log(P(o;lambda)). Choice==1: let sigma_i(beta[t][i])=1 for every t; choice==0: do not scale**/
double backwardstep(double** beta, double a[][N], double b[][L][M], double pi[], int** o, int n, int numt, int choice)
{
    //double scale[T],sum;
    double sum;
    int i,j,t,l,i0;
    //double x[N][T];
    double* scale = (double*)calloc(T, sizeof(double));
    double** x = (double**)malloc(N * sizeof(double*));
    for (int i0 = 0; i0 < N; i0++)
    {
        x[i0] = (double*)calloc(T, sizeof(double));
    }
    for(i=0; i<n; i++)
        for(t=0; t<numt; t++)
            for( l=0, x[i][t]=1; l<L; l++)
            {
                x[i][t]*=b[i][l][o[t][l]];
            }
    for(i=0; i<n; i++)
        beta[numt-1][i]=1;
    if(choice==1)
        scale[numt-1]=scaler(beta[numt-1],n);
    for(t=numt-2; t>=0; t--)
    {
        for(i=0; i<n; i++)
        {
            for(j=0,sum=0; j<n; j++)
                sum+=a[i][j]*x[j][t+1]*beta[t+1][j];
            beta[t][i]=sum;
        }
        if(choice==1)
            scale[t]=scaler(beta[t],n);
    }
    for(i=0,sum=0; i<n; i++)
        sum+=pi[i]*x[i][0]*beta[0][i];
    sum=log(sum);
    if(choice==1)
        for(t=0; t<numt; t++)
            sum+=log(scale[t]);
    free(scale);
    for ( i0 = 0; i0 < N; i0++)
        free(x[i0]);
    free(x);
    return(sum);
}

