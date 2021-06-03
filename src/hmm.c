#include "hmm.h"

//int hmm(int n, int m, int numc, int numt[C], int o[C][T][L],
//        double a[N][N], double b[N][L][M], double pi[N],
//        double alpha[C][T][N], double beta[C][T][N], double gamma[C][T][N],
//        double xi[C][T-1][N][N], double pi_numerator[C][N], double pi_denominator[C],
//        double a_numerator[C][N][N], double a_denominator[C][N],
//        double b_numerator[C][N][L][M], double b_denominator[C][N],
//        double numerator, double denominator,
//        int sigmat, int par, int tn, unsigned int seed,
//        int choicegetobs, int choiceinitial, int choicefb, clock_t timearr[],
//        int *piter, double *plo, double *pl,
//        int maxiter, double eps, int endwithouteps)
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
        int maxiter, double eps, int endwithouteps)

{
    int c,i,j,k,l,iter;
    double *pm[N];
    double logprobold,logprob;
    clock_t t0=clock();


    /**Get observation**/
    //srand(seed);
    //choicegetobs=1;//which kind of observation to be generated
    // genobslen(numt,numc,sigmat,choicegetobs);
    //srand(seed);
    /** for(c=0; c<numc; c++)
         for(int l=0; l<L; l++)
             for(t=0; t<numt[c]; t++) //new
                 o[c][t][l]=rand()%m;  // o[c][t]=rand()%m;//generate observations randomly
    **/

    /**Initialization
    srand(99);
    //choiceinitial=1;//which kind of initialization to be used
    initializer(pi,n,choiceinitial);
    for(i=0; i<n; i++)
    {
        initializer(a[i],n,choiceinitial);
        for(int l=0; l<L; l++)
            initializer(b[i][l],m,choiceinitial);
    }


    /**Print initialization
    for(c=0; c<numc; c++)
    {
        printf("Obersvation %d\n",c);
        for(t=0; t<numt[c]; t++)
        {

            for(l=0; l<L-1; l++)
            {
                printf("%d ",o[c][t][l]);
            }
            printf("%d\n",o[c][t][l]);
        }
        printf("\n");
    }**/
    /** srand(12345);
     //choiceinitial=1;//which kind of initialization to be used
     initializer(pi,n,choiceinitial);
     for(i=0; i<n; i++)
     {
         initializer(a[i],n,choiceinitial);
         for(int l=0; l<L; l++)
             initializer(b[i][l],m,choiceinitial);
     }**/
    /**get pi**/
    FILE *fp2;
    fp2=fopen("E:\\data\\ABpi\\pi.txt","rt");
    if((fp2=fopen("E:\\data\\ABpi\\pi.txt","rt"))==NULL)
    {
        printf("cannot open file pi\n");
        return 0;
    }
    for(i=0; i<n; i++)
    {
        fscanf(fp2,"%lf",&pi[i]);
    }
    fscanf(fp2,"\n");
    fclose(fp2);
    printf("pi:\n");
    for(i=0; i<n; i++)
    {
        printf("%f\t",pi[i]);
    }
    printf("\n");
    /**get A **/
    FILE *fp3;
    fp3=fopen("E:\\data\\ABpi\\A.txt","rt");
    if((fp3=fopen("E:\\data\\ABpi\\A.txt","rt"))==NULL)
    {
        printf("cannot open file A\n");
        return 0;
    }
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
            fscanf(fp3,"%lf",&a[i][j]);
        fscanf(fp3,"\n");
    }
    fclose(fp3);
    printf("A:\n");
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            printf("%f ",a[i][j]);
        }
        printf("\n");
    }
    /**get B **/

    FILE *fp4;
    fp4=fopen("E:\\data\\ABpi\\B.txt","rt");
    if((fp4=fopen("E:\\data\\ABpi\\B.txt","rt"))==NULL)
    {
        printf("cannot open file B\n");
        return 0;
    }
    for(i=0; i<n; i++)
    {
        for(l=0; l<L; l++)
            fscanf(fp4,"%lf ",&b[i][l][1]);
        fscanf(fp4,"\n");
    }
    fclose(fp4);
    printf("B matrix:\n");
    for(i=0; i<n; i++)
        for(int l=0; l<L; l++)
        {
            if(l<L-1)
                printf("%f\t",b[i][l][1]);
            if(l==L-1)
                printf("%f\n",b[i][l][1]);
        }
    printf("\n");
    printf("Initial state distribution pi:\n");
    printvector(pi,n);
    printf("Transision matrix A:\n");
    for(i=0; i<n; i++)
        pm[i]=a[i];
    printmatrix(&pm[0],n,n);//equals----printmatrix(pm,n,n);
    for(i=0; i<n; i++)
        for(int l=0; l<L; l++)
        {
            b[i][l][0]=1-b[i][l][1];
        }
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
//printf("Emission matrix B:\n");
//for(i=0;i<n;i++)
//  pm[i]=b[i];
//printmatrix(pm,n,m);
    printf("Number of observations: %d\n",numc);
    printf("Length of observations:\n");
    for(c=0; c<numc-1; c++)
        printf("%d ",numt[c]);
    printf("%d\n",numt[numc-1]);
    for(c=0; c<numc; c++)
    {
       // for(int l=0; l<L; l++)
       // {
            forwardstep(alpha[c],a,b,pi,o[c],n,numt[c],choicefb);
            backwardstep(beta[c],a,b,pi,o[c],n,numt[c],choicefb);
            getgamma(gamma[c],alpha[c],beta[c],n,numt[c]);
            // getxi(xi[c],alpha[c],beta[c],a,b,o[c],n,numt[c]);
            //  getnumeratordenominator(pi_numerator[c],&pi_denominator[c],a_numerator[c],a_denominator[c],b_numerator[c],b_denominator[c],o[c],gamma[c],xi[c],n,m,numt[c]);
            /** if(iter==500)
             {
                 printf("Number of cores =%d, thread ID =%d, c: %d\n", omp_get_num_threads(), omp_get_thread_num(), c);
             }
            **/
       // }
    }
    return(0);
}


  /**Baum-Welch

    //get logprob
    //choicefb=1;//scale [1] or not [0]
    logprob=0;//omp do not allow for(logprob=0,c=0;c<numc;c++) !!!
    //#pragma omp parallel for reduction(+:logprob)

    //for(c=0; c<numc; c++)
    logprob+=forwardstep(alpha[c],a,b,pi,o[c],n,numt[c],choicefb);
    printf("logprob: %f\n\n",logprob);

    iter=0;
    clock_t tempserial=0;
    clock_t tempparrellel=0;
    clock_t t1=clock();
    do
    {
        clock_t t0=clock();
        //with this statement, the t0 before "do" will be ignored;
       // without this statement, t0 will be the t0 before "do"
        logprobold=logprob;
        //E step
        clock_t t1=clock();
        //#pragma omp parallel for num_threads(tn) schedule(dynamic,1) if(par)
        //for(c=0; c<numc; c++)
        //{
        for(int l=0; l<L; l++)
        {
            //backwardstep(beta[c],a,b,pi,o[c],n,numt[c],choicefb);
            //getgamma(gamma[c],alpha[c],beta[c],n,numt[c]);
            //getxi(xi[c],alpha[c],beta[c],a,b,o[c],n,numt[c]);
            //getnumeratordenominator(pi_numerator[c],&pi_denominator[c],a_numerator[c],a_denominator[c],b_numerator[c],b_denominator[c],o[c],gamma[c],xi[c],n,m,numt[c]);
            backwardstep(beta,a,b,pi,o,n,numt,choicefb);
            getgamma(gamma,alpha,beta,n,numt);
            getxi(xi,alpha,beta,a,b,o,n,numt);
            getnumeratordenominator(pi_numerator,&pi_denominator,a_numerator,a_denominator,b_numerator,b_denominator,o,gamma,xi,n,m,numt);
            //if(iter==500)
             //{
               //  printf("Number of cores =%d, thread ID =%d, c: %d\n", omp_get_num_threads(), omp_get_thread_num(), c);
            // }

        }
        // }
        clock_t t2=clock();
        //if(iter==500)
        //printf("\n");
        //M step
        //update pi
        denominator=0;
        //#pragma omp parallel for reduction(+:denominator) num_threads(tn) schedule(dynamic,1)
        // for(c=0; c<numc; c++)
        denominator+=pi_denominator[c];
        //#pragma omp parallel for???
        for(i=0; i<n; i++)
        {
            numerator=0;
            //#pragma omp parallel for reduction(+:numerator) num_threads(tn) schedule(dynamic,1)
            // for(c=0; c<numc; c++)
            numerator+=pi_numerator[c][i];
            pi[i]=numerator/denominator;
        }
        //update a
        //#pragma omp parallel for???
        for(i=0; i<n; i++)
        {
            denominator=0;
            //for(c=0; c<numc; c++)
            denominator+=a_denominator[c][i];
            for(j=0; j<n; j++)
            {
                numerator=0;
                //for(c=0; c<numc; c++)
                numerator+=a_numerator[c][i][j];
                a[i][j]=numerator/denominator;
            }
        }
        //update b
        //#pragma omp parallel for???
        for(i=0; i<n; i++)
        {
            denominator=0;
            //for(c=0; c<numc; c++)
            denominator+=b_denominator[c][i];
            for(int l=0; l<L; l++)
            {
                for(k=0; k<m; k++)
                {
                    numerator=0;
                    //  for(c=0; c<numc; c++)
                    numerator+=b_numerator[c][i][l][k];
                    b[i][l][k]=numerator/denominator;
                }
            }
        }
        //get logprob
        logprob=0;
        clock_t t3=clock();
        //#pragma omp parallel for reduction(+:logprob) num_threads(tn) schedule(dynamic,1) if(par)
        //for(c=0; c<numc; c++)
        //{
        logprob+=forwardstep(alpha[c],a,b,pi,o[c],n,numt[c],choicefb);
        //if(iter==500)
       // {
          //  printf("Number of cores =%d, thread ID =%d, c: %d\n", omp_get_num_threads(), omp_get_thread_num(), c);
      //  }
        // }
        clock_t t4=clock();
        iter++;
        //Print iterations
        //if(iter%1000==0||iter<5)
        if(iter<300)
            printf("Iteration: %d\n",iter);
        //if(iter%1000==0||iter<5)
        if(iter<300)
        {
            printf("Initial state distribution pi:\n");
            printvector(pi,n);
            printf("Transision matrix A:\n");
            for(i=0; i<n; i++)
                pm[i]=a[i];
            printmatrix(pm,n,n);
            printf("Emission matrix B:\n");
            for(int l=0; l<L; l++)
            {
                printf("The %d of B matrix:\n",l);
                for(i=0; i<n; i++)
                    for(int k=0; k<m; k++)
                    {
                        if(k<m-1)
                            printf("%f\t",b[i][l][k]);
                        if(k==m-1)
                            printf("%f\n",b[i][l][k]);
                    }
                printf("\n");
            }
            printf("logprob: %.10f\n",logprob);
            printf("logprob - old_logprob: %.10f\n\n",logprob-logprobold);
        }
        clock_t t5=clock();
        tempserial+=t1-t0+t3-t2+t5-t4;
        tempparrellel+=t2-t1+t4-t3;
    }
    while(iter<maxiter&&((logprob-logprobold)>eps||endwithouteps)); //endwithouteps==0, then it equals iter<maxiter&&(logprob-logprobold)>eps; endwithouteps==1, then it equaks iter<maxiter
    clock_t t2=clock();
    timearr[0]=t2-t0;
    timearr[1]=t1-t0;
    timearr[2]=t2-t1;
    timearr[3]=tempserial;
    timearr[4]=tempparrellel;
    *piter=iter;
    *plo=logprobold;
    *pl=logprob;
    return(0);
    **/



