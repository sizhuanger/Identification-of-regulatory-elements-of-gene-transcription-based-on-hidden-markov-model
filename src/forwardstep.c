#include "hmm.h"

/**Return log(P(o;lambda)). Choice==1: let sigma_i(alpha[t][i])=1 for every t; choice==0: do not scale**/
double forwardstep(double** alpha, double a[][N], double b[][L][M], double pi[], int** o, int n, int numt, int choice)
{
    //printf("o[0]:%d", o[0][0]);
    //system("pause");
    //double scale[T],sum;
    double sum;
    //double scale[T];
    double* scale = (double*)calloc(T, sizeof(double));
    int i,j,l,t,i0;
    //double x[N][T];
    double** x = (double**)malloc(N * sizeof(double*));
    for (int i0 = 0; i0 < N; i0++)
    {
        x[i0] = (double*)calloc(T, sizeof(double));
    }
    for(i=0; i<n; i++)
        for(t=0; t<numt; t++)
            for( l=0,x[i][t]=1; l<L; l++)
            {
                //printf("%d", o[t][l]);
                x[i][t]*=b[i][l][o[t][l]];
                //system("pause");
            }
    for (i=0;i<n; i++)
    {
        alpha[0][i] = pi[i]*x[i][0];
    }
    if(choice==1)
        scale[0]=scaler(alpha[0],n);

    for(t=1; t<numt; t++)
    {

        for(j=0; j<n; j++)
        {
            for(i=0,sum=0; i<n; i++)
                sum+=alpha[t-1][i]*a[i][j];
            sum*=x[j][t];
            alpha[t][j]=sum;

        }
        if(choice==1)
            scale[t]=scaler(alpha[t],n);

    }
    if(choice==1)
        for(t=0,sum=0; t<numt; t++)
            sum+=log(scale[t]);
    else
    {
        for(i=0,sum=0; i<n; i++)
            sum+=alpha[numt-1][i];
        sum=log(sum);
    }
    free(scale);
    for ( i0 = 0; i0 < N; i0++)
    {
        free(x[i0]);
    }
    free(x);
    return(sum);
}

