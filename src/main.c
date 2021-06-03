#include "hmm.h"

int main()
{
    int m;//number of possible values of observations
    int n;//number of states
    int numc;//number of observed sequences
    static int numt[C];//length observed sequences
    //static int o[C][T][L];//observed sequences
    static double a[N][N];//transision matrix
    static double b[N][L][M];//emission matrix
    static double pi[N];//initial state distribution
    //static double alpha[C][T][N];//alpha[t][i]=P(o_0,o_1,...,o_t,s_t=i;lambda),t=0,...,T-1,i=0,...,n-1
    //static double beta[C][T][N];//beta[t][i]=P(o_t+1,...,o_T-1|s_t=i;lambda),t=0,...,T-1,i=0,...,n-1
    //static double gamma[C][T][N];//gamma[t][i]=P(s_t=i|o;lambda),t=0,...,T-1,i=0,...,n-1
    //static double xi[C][T-1][N][N];//xi[t][i][j]=P(s_t=i,s_t+1=j|o,lambda),t=0,...,T-2,i=0,...,n-1,j=0,...,n-1
    static double pi_numerator[C][N];
    static double pi_denominator[C];
    static double a_numerator[C][N][N];
    static double a_denominator[C][N];
    static double b_numerator[C][N][L][M];
    static double b_denominator[C][N];
    static double numerator;
    static double denominator;
    static double *pm[N];//used for printmatrix()
    //////////////////////////////////////////////////////////////////
    //static int o[C][T][L];//observed sequences
    //static double alpha[C][T][N];//alpha[t][i]=P(o_0,o_1,...,o_t,s_t=i;lambda),t=0,...,T-1,i=0,...,n-1
    //static double beta[C][T][N];//beta[t][i]=P(o_t+1,...,o_T-1|s_t=i;lambda),t=0,...,T-1,i=0,...,n-1
    //static double gamma[C][T][N];//gamma[t][i]=P(s_t=i|o;lambda),t=0,...,T-1,i=0,...,n-1
    //static double xi[C][T - 1][N][N];//xi[t][i][j]=P(s_t=i,s_t+1=j|o,lambda),t=0,...,T-2,i=0,...,n-1,j=0,...,n-1
    //int hiddeni[T];
    //static int o[C][T][L];
    int*** o = (int***)malloc(C * sizeof(int**));
    for (int i0 = 0; i0 < C; i0++)
    {
        o[i0] = (int**)malloc(T * sizeof(int*));
        for (int i1 = 0; i1 < T; i1++)
        {
            o[i0][i1] = (int*)calloc(L,sizeof(int));
        }
    }
    for (int i0 = 0; i0 < C; i0++)
    {
        for (int i1 = 0; i1 < T; i1++)
        {
            for (int i2 = 0; i2 < L; i2++)
            {
                //printf("%d", o[i0][i1][i2]);
            }
        }
    }

    /*static double alpha[C][T][N];*/
    double*** alpha=(double***)malloc(C*sizeof(double**));
    for (int i0 = 0; i0 < C; i0++)
    {
        alpha[i0] = (double**)malloc(T * sizeof(double*));
        for (int i1 = 0; i1 < T; i1++)
        {
            alpha[i0][i1] = (double*)calloc(N, sizeof(double));
        }
    }


    //static double beta[C][T][N];
    double*** beta = (double***)malloc(C * sizeof(double**));
    for (int i0 = 0; i0 < C; i0++)
    {
        beta[i0] = (double**)malloc(T * sizeof(double*));
        for (int i1 = 0; i1 < T; i1++)
        {
            beta[i0][i1] = (double*)calloc(N, sizeof(double));
        }
    }
    //static double gamma[C][T][N];
    double*** gamma = (double***)malloc(C * sizeof(double**));
    for (int i0 = 0; i0 < C; i0++)
    {
        gamma[i0] = (double**)malloc(T * sizeof(double*));
        for (int i1 = 0; i1 < T; i1++)
        {
            gamma[i0][i1] = (double*)calloc(N, sizeof(double));
        }
    }
    //static double xi[C][T - 1][N][N];
    double**** xi = (double****)malloc(C * sizeof(double***));
    for (int i0 = 0; i0 < C; i0++)
    {
        xi[i0] = (double***)malloc((T - 1) * sizeof(double**));
        for (int i1 = 0; i1 < T - 1; i1++)
        {
            xi[i0][i1] = (double**)malloc(N * sizeof(double*));
            for (int i2 = 0; i2 < N; i2++)
            {
                xi[i0][i1][i2] = (double*)calloc(N, sizeof(double));
            }
        }
    }
    //int hiddeni[T]
    int* hiddeni = (int*)calloc(T, sizeof(int));

    //////////////////////////////////////////////////////////////////
    unsigned int seed;
    int sigmat,par,tn,choicegetobs,choiceinitial,choicefb,maxiter,endwithouteps,size;
    clock_t timearr[5];
    int i,l,c,itn=0,jsize,t;
    int iter,*piter=&iter;
    double logprobold,logprob,*plo=&logprobold,*pl=&logprob,eps;
    clock_t extensibility[11][10];

    char *file[]= {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"};
    char qianzhui[50]="E:\\data\\H1_chr";
    char houzhui[20]="_binary_test.txt";
    char final[50]="";
    char finally[50]="";
    char out[20]="_out_test.txt";

    maxiter=300;
    eps=1e-3;
    endwithouteps=0;

    //for(itn=0; itn<11; itn++)
    // for(jsize=0; jsize<10; jsize++)
    //{
    jsize=100;
    par=1;
    if(itn==0) par=0;
    tn=itn;//set num_threads
    size=jsize+1;

    /**Given n and m**/
    n=15;
    m=2;
    /**Get observation**/
    numc=1;
    sigmat=200*size;
    choicegetobs=0;//which kind of observation to be generated
    /**Initialization**/
    choiceinitial=1;//which kind of initialization to be used
    /**Scale or not**/
    choicefb=1;//scale [1] or not [0]
    //seed=2;

    for(c=0; c<numc; c++)
    {
        strcpy(final,"");
        strcat(final,qianzhui);
        strcat(final,file[c]);
        strcat(final,houzhui);
        FILE *fp;
        fp=fopen(final,"rt");
        if((fp=fopen(final,"rt"))==NULL)
        {
            printf("cannot open file\n");
            return 0;
        }
        for (int x=0; x<2; x++)
            fscanf(fp,"%*[^\n]%*c");
        t=0;
        while (!feof(fp))
        {

            for(l=0; l<L; l++)
            {
                fscanf(fp,"%d",&o[c][t][l]);
            }
            fscanf(fp,"\n");
            t++;
        }
        numt[c]=t;
        fclose(fp);
        printf("Observation%d:\n",c+1);


        /** for(j=0; j<numt[c]; j++)
         {
             for(l=0; l<L; l++)
                 printf("%d\t",o[c][j][l]);
             printf("\n");
         }**/
    }
    /**Train hmm**/
    hmm(n, m, numc, numt, o,
        a, b, pi,
        alpha, beta, gamma,
        xi, pi_numerator, pi_denominator,
        a_numerator, a_denominator,
        b_numerator, b_denominator,
        numerator, denominator,
        sigmat, par, tn, seed,
        choicegetobs, choiceinitial, choicefb, timearr,
        piter, plo, pl,
        maxiter, eps, endwithouteps);

    /**Print the result**/
    printf("-------------------------------Result-------------------------------\n");
    //printf("Observations:\n");
    //for(c=0;c<numc;c++)
    //    printobservation(o[c],numt[c]);
    for(c=0; c<numc; c++)
    {
        // if(c<numc-1)
        //{
        printf("chr%d hiddenstates:\n",c+1);
        gethiddenstates(hiddeni,gamma[c], numt[c],n);
        strcpy(finally,"");
        strcat(finally,qianzhui);
        strcat(finally,file[c]);
        strcat(finally,out);
        FILE *fp1;
        fp1=fopen(finally,"wt");
        if(NULL==fp1)
        {
            printf("open error\n");
        }
        else
        {
            for(t=0; t<numt[c]; t++)
            {
                fprintf(fp1,"%d",hiddeni[t]);
                fprintf(fp1,"\n");
            }
        }
        // printf("Iteration: %d\n",iter);
        printf("Number of states: %d\n",n);
        printf("Number of possible values of observations: %d\n",m);
        printf("Initial state distribution pi:\n");
        printvector(pi,n);
        printf("Transision matrix A:\n");
        for(i=0; i<n; i++)
            pm[i]=a[i];
        printmatrix(&pm[0],n,n);//equals----printmatrix(pm,n,n);
        printf("Emission matrix B:\n");
        printf("when observation is 1 the B:\n");
        for(i=0; i<n; i++)
            for(int l=0; l<L; l++)
            {

                if(l<L-1)
                    printf("%f\t",b[i][l][1]);
                if(l==L-1)
                    printf("%f\n",b[i][l][1]);
            }
        printf("when observation is 0 the B:\n");
        for(i=0; i<n; i++)
            for(int l=0; l<L; l++)
            {

                if(l<L-1)
                    printf("%f\t",b[i][l][0]);
                if(l==L-1)
                    printf("%f\n",b[i][l][0]);
            }
        printf("\n");
        printf("Number of observed sequences: %d\n",numc);
        printf("Length of observations:\n");
        for(c=0; c<numc-1; c++)
            printf("%d ",numt[c]);
        printf("%d\n",numt[numc-1]);
        printf("logprob: %.10f\n",logprob);
        printf("logprob - old_logprob: %.10f\n\n",logprob-logprobold);

        /**Time used**/

        printf("Total time = %d ms\n",timearr[0]);
        printf("|--Initialization time = %d ms\n",timearr[1]);
        printf("|--Iteration time = %d ms\n",timearr[2]);//show the time used
        printf("   |--Serial time in iterations = %d ms\n",timearr[3]);
        printf("   |--Parallel time in iterations = %d ms\n",timearr[4]);

        extensibility[itn][jsize]=timearr[2];
        //printf("%d\n",extensibility[itn][jsize]);
        //getchar();
        // }
        for(itn=0; itn<11; itn++)
        {
            for(jsize=0; jsize<10-1; jsize++)
            {
                printf("%d\t",extensibility[itn][jsize]);
            }
            printf("%d\n",extensibility[itn][jsize]);
        }
        return 0;
    }
    }
