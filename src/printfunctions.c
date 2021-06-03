#include "hmm.h"

int printmatrix(double **mat, int nrow, int ncol)
{
    int i;
    int printvector(double arr[], int len);
    for(i=0; i<nrow; i++)
        printvector(mat[i],ncol);
    return(0);
}

int printvector(double arr[], int len)
{
    int i;
    for(i=0; i<len-1; i++)
        printf("%f\t",arr[i]);
    printf("%f\n",arr[i]);
    return(0);
}


