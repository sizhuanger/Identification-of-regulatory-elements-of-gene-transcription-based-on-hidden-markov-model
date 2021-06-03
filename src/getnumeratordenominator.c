#include "hmm.h"

int getnumeratordenominator(double pi_numerator[], double *pi_denominator, double a_numerator[][N],
                            double a_denominator[], double b_numerator [][L][M], double b_denominator[],
                            int** o, double** gamma, double*** xi,int n,int m,int numt)
{
    int i,j,k,t;
    for(i=0; i<n; i++)
        pi_numerator[i]=gamma[0][i];
    *pi_denominator=1;
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
            for(t=0,a_numerator[i][j]=0; t<numt-1; t++)
                a_numerator[i][j]+=xi[t][i][j];
        for(t=0,a_denominator[i]=0; t<numt-1; t++)
            a_denominator[i]+=gamma[t][i];
        b_denominator[i]=a_denominator[i]+gamma[numt-1][i];
        for(int l=0; l<L; l++)
        {
            for(k=0; k<m; k++)
                for(t=0,b_numerator[i][l][k]=0; t<numt; t++)
                    if(o[t][l]==k)
                        b_numerator[i][l][k]+=gamma[t][i];
        }
    }
    return(0);
}
