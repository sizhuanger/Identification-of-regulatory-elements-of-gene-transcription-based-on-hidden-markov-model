#include "hmm.h"

//int getgamma(double gamma[][N], double alpha[][N], double beta[][N], int n, int numt)
int getgamma(double** gamma, double** alpha, double** beta, int n, int numt)
{
	//for (int i0 = 0; i0 < 10; i0++) {
	//	printf("tt:%lf", beta[i0][i0]);
	//}
    int i,t;
    for(t=0; t<numt; t++)
    {
        for(i=0; i<n; i++)
            gamma[t][i]=alpha[t][i]*beta[t][i];
        scaler(gamma[t],n);
    }
    return(0);
}

//int gethiddenstates(int hiddeni[],double gamma[][N],int numt,int n)
int gethiddenstates(int hiddeni[],double **gamma,int numt,int n)
{
    int i,t;
    double maxgamma;
    for(t=0; t<numt; t++)
    {
        maxgamma=0;
        for(i=0; i<n; i++)
        {

            if(gamma[t][i]>maxgamma)
            {
                maxgamma=gamma[t][i];
                hiddeni[t]=i;
            }
        }
    }
    return(0);
}
