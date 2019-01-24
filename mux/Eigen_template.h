/*Intel® Math Kernel Library LAPACK Examples
/*******************************************************************************
*  Copyright (C) 2009-2015 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
* https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/ssyevr_ex.c.htm
********************************************************************************
*/

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<random>
#include<ctime>
//#include<algorithm>
//#include<iterator>
#include<math.h>
#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>
#include<chrono>
#include<ctime>
#include<cmath>


//Paper: Sarika Jalan and Priodyuti Pradhan, Localization of multilayer networks by optimized single-layer rewiring, Phys. Rev. E 97, 042314 (2018)
// Written: Priodyuti Pradhan, Complex Systems Lab, IIT Indore
// contact: priodyutipradhan@gmail.com  

using namespace std::chrono;

extern "C"{ void ssyevr_(char* jobz, char* range, char* uplo, int* n, float* a,
                       int* lda, float* vl, float* vu, int* il, int* iu, float* abstol,
                       int* m, float* w, float* z, int* ldz, int* isuppz, float* work,
                       int* lwork, int* iwork, int* liwork, int* info );
          };


template<class X>double Eigen11(X **g, int N, double Ev[]);
template<class X>void print_matrix(int m, int n, X* a, int lda );
//double IPR_MUX(double *Ev, int N, double *ipr);
//double ipr_layer_contribution_mux(double *Ev, unsigned int N, unsigned int l, double *ipr);
double IPR(double Ev[],int N);


double IPR(double *Ev,int N)
{
 double Num = 0.0;
 for(int i=0;i<N;i++)
    Num = Num + pow(Ev[i],4);   
 return Num;
}

/*double IPR_MUX(double *Ev,int N,double *ipr)
{
 double Num = 0.0;
 for(int i=0;i<N/2;i++)
    Num = Num + pow(Ev[i],4);
 ipr[1] = Num;
 Num = 0.0;
 for(int i=N/2;i<N;i++)
    Num = Num + pow(Ev[i],4);
 ipr[2] = Num;     
 
 return 1.0;
}*/

/*double ipr_layer_contribution_mux(double *Ev,unsigned int N,unsigned int l,double *ipr)
{
 double Num = 0.0;
 for(int i=0;i<N/2;i++)
    Num = Num + pow(Ev[i],4);
 ipr[1] = Num;
 Num = 0.0;
 for(int i=N/2;i<N;i++)
    Num = Num + pow(Ev[i],4);
 ipr[2] = Num;     
 
 return 1.0;
*/

template<class X>double Eigen11(X **g,int N, double Ev[])
 {
     int n = N, il, iu, m, lda = N, ldz = N, info, lwork, liwork;
     float* a;
     a = new float[N*N];
     //std::vector<double> v(N*N);
     unsigned long int Val = N*N;
     for(unsigned long int i=0;i<Val;i++)
       a[i] = 0.0;
     unsigned long int k=0;
     int M = 1; 
     for(int i = 0;i<N;i++)
     {
       for(int j = 0;j<M;j++)
       {
         a[k++] = g[i][j];
       }
       k = k+N-M;
       M++;
     }

      /*for(unsigned long int i=0;i<Val;i++)
        cout<<a[i]<<' ';
      cout<<endl; */
        float abstol, vl, vu;
        int iwkopt;
        int* iwork;
        float wkopt;
        float* work;
        /* Local arrays */
        int* isuppz;
        isuppz = new int[N];    
        float* z;
        z = new float[N];
          
        float w[1];   
        
        /* Executable statements */
        //cout<<" DSYEVR Example Program Results"<<endl;
        /* Negative abstol means using the default value */
        abstol = -1.0;
        /* Set il, iu to compute NSELECT smallest eigenvalues */
        il = N;
        iu = N;
        /* Query and allocate the optimal workspace */
        lwork = -1;
        liwork = -1;
        char JOBZ = 'V';
        char RANGE = 'I';
        char UPLO = 'U';
                
        ssyevr_(&JOBZ, &RANGE, &UPLO, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, 
                                                                        &wkopt, &lwork, &iwkopt, &liwork, &info );       
 
        lwork = (int)wkopt;
        work = (float*)malloc( lwork*sizeof(float) );
        liwork = iwkopt;
        iwork = (int*)malloc( liwork*sizeof(int) );
        /* Solve eigenproblem */
                
        ssyevr_(&JOBZ, &RANGE, &UPLO, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, 
                                                                           work, &lwork, iwork, &liwork, &info );
        /* Check for convergence */
        if( info > 0 ) {
               cout<< "The algorithm failed to compute eigenvalues."<<endl;
        }
        for(unsigned int i=0;i<N;i++)
            Ev[i] = z[i];
        delete [] a;
        delete [] z;
        delete [] isuppz;
       
        free( (void*)iwork );
        free( (void*)work );
        return w[0];

} /* End of SSYEVR Example */


