#include "hmm.h"

/**choice==0: prob[i]=1/n; choice==1: prob[i]=rand()%RANDMAX;**/
int initializer(double prob[], int n, int choice)
{
    int i,sum;
    if(choice==0)
    {
        for(i=0; i<n; i++)
            prob[i]=1.0/n;
        return(0);
    }
    if(choice==1)
    {
        for(i=0,sum=0; i<n; i++)
        {
            prob[i]=rand()%RANDMAX+1;
            sum+=prob[i];
        }
        for(i=0; i<n; i++)
            prob[i]/=(double)sum;
        return(1);
    }
    return(-1);//choice!=0||1, have not done iitialization
}

int genobslen(int numt[], int numc, int totallen, int choice)
{
    int c,temp;
    if(choice==0)//generate length of observations randomly
    {
        temp=totallen/numc;
        for(c=0; c<numc; c++)
            numt[c]=rand()%(temp-1)+2;
        return(0);
    }
    if(choice==1)//generate balanced length of observations
    {
        balance(numt,0,numc,totallen);
        return(1);
    }
    if(choice==2)//generate unbalanced length of observations
    {
        numt[0]=totallen/2;
        balance(numt,1,numc,totallen-totallen/2);
        return(2);
    }
    return(-1);
}

int balance(int arr[], int first, int last, int total)
{
    int i,remainder,fold,len;
    len=last-first;
    remainder=total%len;
    fold=total/len;
    for(i=first; i<first+remainder; i++)
        arr[i]=fold+1;
    /**The following sentence trouble me a lot! Be careful about array bound!**/
    for(; i<last; i++) //i<last cannot be i<=last, or it will change numt[C] (out of bound!), which can be a global variant in memory, see numc!!!
        arr[i]=fold;
    return(0);
}
