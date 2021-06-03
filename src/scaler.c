#include "hmm.h"

/**Let sum equal one and return the sum**/
double scaler(double *arr, int len)
{
    double sum;
    int i;
    for(i=0,sum=0; i<len; i++)
        sum+=arr[i];
    for(i=0; i<len; i++)
        arr[i]/=sum;
    return(sum);
}
