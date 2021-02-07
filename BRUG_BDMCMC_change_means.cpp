/***************************************************************************************
*     Copyright (C) 2020 Chan Ga Ming Angus
*
*     Maintainer: Chan Ga Ming Angus <chan.ga.ming.angus@gmail.com>
*     Source: https://github.com/Anguscgm/BRUG_BDMCMC
***************************************************************************************/

// Compile with:
// g++ BRUG_BDMCMC_change_means.cpp -llapack -lblas -std=c++11 -o BRUG_BDMCMC_change_means

/***************************************************************************************
*    Reference
****************************************************************************************
*    Title: BDgraph source code
*    Author: Reza Mohammadi, Ernst C. Wit
*    Date: 2012 - 2020
*    Code version: 2.63
*    Availability: https://github.com/cran/BDgraph
*
***************************************************************************************/

#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <cstring>
#include <getopt.h>
#include <algorithm>
using namespace std;

extern "C" void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);
extern "C" void dtrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);
extern "C" void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);
extern "C" void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double* C, int* LDC);
extern "C" void dposv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, int* INFO);
extern "C" void dsymv_(char* UPLO, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
extern "C" void dsyr_(char* UPLO, int* N, double* ALPHA, double* X, int* INCX, double* A, int* LDA);
extern "C" double ddot_(int* N, double* DX, int* INCX, double* DY, int* INCY);
extern "C" void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);

void update_alpha(double alpha0[], double alpha1[], double K[], double X[], double Y_bar[], double mu[], int ns[], int p, int L)
{
    double mu_alpha0 = 0, mu_alpha1 = 0;
    double sigma2_alpha0 = 1; double sigma2_alpha1 = 1;
    double sigma2_MH_alpha0 = 0.1; double sigma2_MH_alpha1 = 0.01;
    default_random_engine generator;
    double alpha_old, alpha_new, logr;
    char transN = 'N';
    int one = 1; double d_one = 1.0, d_zero = 0.0;
    vector<double> S_diff(p*p);
    vector<double> SK(p*p);
    double trace;
    for (int i=0; i<p; i++)
    {
        int updated = 0;
        fill(&S_diff[0], &S_diff[p*p], 0);
        // Update alpha0
        alpha_old = alpha0[i];
        normal_distribution<double> N(alpha_old,sigma2_MH_alpha0); //Can change the sd later.
        alpha_new = N(generator);
        logr = (pow(alpha_old-mu_alpha0,2.0)-pow(alpha_new-mu_alpha0,2.0))/2/sigma2_alpha0;
        for (int l=0; l<L; l++)
        {
            S_diff[i*(p+1)] = pow(alpha_new,2.0)-pow(alpha_old,2.0) - 2*(alpha_new-alpha_old)*(Y_bar[l*p+i]-alpha1[i]*X[l]);
            for (int j=0; j<p && j!=i; j++)
            {
                S_diff[i*p+j] = (alpha_old-alpha_new)*(Y_bar[l*p+j]-mu[l*p+j]);
                S_diff[j*p+i] = S_diff[i*p+j];
            }
            dgemm_( &transN, &transN, &p, &p, &p, &d_one, &S_diff[0], &p, &K[l*p*p], &p, &d_zero, &SK[0], &p );
            trace = 0; //reset
            for (int j=0; j<p; j++) trace += SK[j*(p+1)];
            logr -= ns[l]*trace/2.0;
        }
        if (logr > log((double)rand()/RAND_MAX))
        {
            updated = 1;
            alpha0[i] = alpha_new;
        }
        
        // alpha1
        alpha_old = alpha1[i];
        normal_distribution<double> N_2(alpha_old,sigma2_MH_alpha1); //Can change the sd later.
        logr = (pow(alpha_old-mu_alpha1,2.0)-pow(alpha_new-mu_alpha1,2.0))/2/sigma2_alpha1;
        for (int l=0; l<L; l++)
        {
            S_diff[i*(p+1)] = (pow(alpha_new,2.0)-pow(alpha_old,2.0))*pow(X[l],2.0) -2*(alpha_new-alpha_old)*(Y_bar[l*p+i]-alpha0[i]);
            for (int j=0; j<p && j!=i; j++)
            {
                S_diff[i*p+j] = (alpha_old-alpha_new)*X[l]*(Y_bar[l*p+j]-mu[l*p+j]);
                S_diff[j*p+i] = S_diff[i*p+j];
            }
            dgemm_( &transN, &transN, &p, &p, &p, &d_one, &S_diff[0], &p, &K[l*p*p], &p, &d_zero, &SK[0], &p );
            trace = 0; //reset
            for (int j=0; j<p; j++) trace += SK[j*(p+1)];
            logr -= ns[l]*trace/2.0;
        }
        if (logr > log((double)rand()/RAND_MAX))
        {
            updated = 1;
            alpha1[i] = alpha_new;
        }
        if (updated)
        {
            for (int l = 0; l<L; l++)
                mu[l*p + i] = alpha0[i] + alpha1[i]*X[l];
        }
    }
}

// Write a function to update \beta_lm
void update_beta(double beta0[], double beta1[], int G[], double X[], int p, int L, int ns[])
{
    int lm = 0;
    default_random_engine generator;
    for (int i=1; i<p; i++) // all potential edges.
        for (int j = 0; j<i; j++)
        {
            int ij = i*p + j;
            // Update beta0
            //prior N(-1000, 1.0)
            //MH N(beta^{t-1}, 0.1)
            double beta_old = beta0[lm];
            normal_distribution<double> N(beta_old,0.01); //Can change the sd later.
            double beta_new = N(generator);
            double logr = (pow(beta_old + 1000, 2.0) - pow(beta_new + 1000, 2.0))/2; //assuming prior of -1000 for beta0
            for (int l=0; l<L; l++) //groups
            {
                double x = X[l];
                logr += ns[l]*(G[l*p*p + ij]*(beta_new - beta_old) - log(1 + exp(beta_new + beta1[lm]*x)) + log(1 + exp(beta_old + beta1[lm]*x)));
            }
            if (logr > log((double)rand()/RAND_MAX)) beta0[lm] = beta_new;
            
            // Update beta1
            //prior N(0, 1.0)
            //MH N(beta^{t-1}, 0.1)
            beta_old = beta1[lm];
            normal_distribution<double> N_2(beta_old,0.01); //Can change the sd later.
            beta_new = N_2(generator);
            //prior
            logr = (pow(beta_old, 2.0) - pow(beta_new, 2.0))/2; // assume \sigma  = 1.0?
            for (int l=0; l<L; l++) //groups
            {
                double x = X[l];
                logr += ns[l]*(G[l*p*p + ij]*(beta_new - beta_old)*x - log(1 + exp(beta0[lm] + beta_new*x)) + log(1 + exp(beta0[lm] + beta_new*x)));
            }
            if (logr > log((double)rand()/RAND_MAX)) beta1[lm] = beta_new;
            lm++;
        }
}

int factorial(int n)
{
    int result = 1;
    for( int i=2; i<=n; i++) result *= i;
    return(result);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//  Determinant of a symmetric possitive-definite matrix ( A )
//       > > > > > > > > >  WARNING: Matrix you pass is overwritten < < < < < < < < < 
//  For any symmetric PD Matrix A, we have: |A| = |T| ^ 2, where T is cholesky decomposition of A. 
//  Thus, |T| = \prod_{i = 1}^p T_{ii}.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
double determinant( double A[], int dim )
{
	char uplo = 'U';
	int info, dim1 = dim + 1;
	
	dpotrf_( &uplo, &dim, &A[0], &dim, &info );

	double result = 1;
	for( int i = 0; i < dim; i++ ) result *= A[ i * dim1 ];
	
	return ( result * result );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Part of function "gnorm"
// which is for calculating Normalizing constant of G-Wishart distribution 
// based on Monto Carlo algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_exp_mc( int G[], int nu[], int b_c, double H[], int check_H, int mc_iter, int dim, double f_T[] )
{
	int iter, i, j, ij, h, r, pxp = dim * dim;
	
	double sumPsi, sumPsiH, sumPsiHi, sumPsiHj;
	double max_numeric_limits_ld = numeric_limits<double>::max() / 1000;
	double min_numeric_limits_ld = numeric_limits<double>::min() * 1000;
	
	vector<double> psi( pxp, 0.0 );      

	//GetRNGstate();
    default_random_engine generator;
    normal_distribution<double> rnorm( 0.0, 1.0 );
	if( check_H == 1 )
	{ 
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
            {
                gamma_distribution<double> rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 );
				psi[ i * dim + i ] = sqrt( rgamma(generator));
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );
            }
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); else psi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? rnorm(generator) : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						psi[ ij ] = 0.0; // it's not necessary
						if( i > 0 )  
						{
							sumPsi = 0.0;
							//sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] )
							// for( h = 0; h < ( i - 1 ); h++ )
							for( h = 0; h < i; h++ )
							{
								if( sumPsi > max_numeric_limits_ld ) sumPsi = max_numeric_limits_ld;	
								if( sumPsi > min_numeric_limits_ld ) sumPsi = min_numeric_limits_ld;	
								sumPsi += ( psi[ i * dim + h ] * psi[ j * dim + h ] );
							}
							
							//psi[i, j] <- - sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] ) / psi[i, i]
							psi[ ij ] = - sumPsi / psi[ i * dim + i ];
						}
						
						if( psi[ ij ] > max_numeric_limits_ld ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] > min_numeric_limits_ld ) psi[ ij ] = min_numeric_limits_ld;	
												
						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] > max_numeric_limits_ld ) f_T[ iter ] = max_numeric_limits_ld;			
		} 
	}else{
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
            {
                gamma_distribution<double> rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 );
				psi[ i * dim + i ] = sqrt(rgamma( generator ));
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );
            }
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); elsepsi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? rnorm(generator) : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						//psi[i, j] = - sum( psi[ i, i : ( j - 1 ) ] * H[ i : ( j - 1 ), j ] )
						sumPsiH = 0.0;
						for( h = i; h < j; h++ )
						{
							if( sumPsiH > max_numeric_limits_ld ) sumPsiH = max_numeric_limits_ld;	
							if( sumPsiH > min_numeric_limits_ld ) sumPsiH = min_numeric_limits_ld;	
							sumPsiH += ( psi[ h * dim + i ] * H[ j * dim + h ] ); 
						}
						psi[ ij ] = - sumPsiH;
						
						if( i > 0 )  //if( i > 1 )
							for( r = 0; r < i; r++ ) //for( r in 1 : ( i - 1 ) )
							{
								//sum( psi[ r, r : i ] * H[ r : i, i ] )
								sumPsiHi = 0.0;
								for( h = r; h < i + 1; h++  )
								{
									if( sumPsiHi > max_numeric_limits_ld ) sumPsiHi = max_numeric_limits_ld;	
									if( sumPsiHi > min_numeric_limits_ld ) sumPsiHi = min_numeric_limits_ld;	
									sumPsiHi += ( psi[ h * dim + r ] * H[ i * dim + h ] );	
								}
									
								//sum( psi[ r, r : j ] * H[ r : j, j ] ) )
								sumPsiHj = 0.0;
								for( h = r; h < j + 1; h++  )
									sumPsiHj += ( psi[ h * dim + r ] * H[ j * dim + h ] );
								
								//psi[i, j] <- psi[i, j] - ( ( sum( psi[ r, r : i ] * H[ r : i, i ] ) ) * ( sum( psi[ r, r : j ] * H[ r : j, j ] ) ) ) / ( psi[i, i] )
								psi[ ij ] -= ( sumPsiHi * sumPsiHj ) / psi[ i * dim + i ];
							}

						if( psi[ ij ] > max_numeric_limits_ld ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] > min_numeric_limits_ld ) psi[ ij ] = min_numeric_limits_ld;	

						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] > max_numeric_limits_ld ) f_T[ iter ] = max_numeric_limits_ld;			
		}
	}
	//PutRNGstate();	
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse( double A[], double A_inv[], int dim )
{
    int info;
    char uplo = 'U';

    // creating an identity matrix
    #pragma omp parallel for
    for( int i = 0; i < dim; i++ )
        for( int j = 0; j < dim; j++ )
            A_inv[ j * dim + i ] = ( i == j );
    
    // LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
    //F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info FCONE );
    dposv_(&uplo, &dim, &dim, A, &dim, A_inv, &dim, &info);
}

double gnorm(int G[], int b, double D[], int p, int iter)
{
    double result;
    
    // get Ti
    vector<double> Ti(p*p);
    char uplo = 'U'; int info;
    inverse(D, &Ti[0], p);
    dpotrf_( &uplo, &p, &Ti[0], &p, &info );
    
    vector<double> H(p*p);
    int check_H = 1;
    for (int i=0; i<p; i++)
        for (int j=0; j<p; j++)
            H[i*p + j] = Ti[i*p + j]/Ti[i*p+i];
    int qp = p*(p-1)/2;
    for (int i=1; i<p; i++)
    {
        for (int j=0; j<i; j++)
        {
            if(H[i*p + j] != 0)
            {
                check_H = 0;
                break;
            }
        }
        if (check_H == 0) break;
    }
    
    vector<int> nu(p, 0);
    for (int i=0; i<p; i++) for (int j=i; j<p; j++) nu[i] += G[i*p + j];
    int size_graph = 0; for (int i=0; i<p; i++) size_graph += nu[i];
    
    if (size_graph == qp)
    {
        double sum_lgamma_bnu = 0;
        for (int i=0; i<p; i++) sum_lgamma_bnu += lgamma((b + nu[i])/(double)2);
        result = ( (double)size_graph / 2 ) * log( M_PI ) + ( (double)p * ( b + p - 1 ) / 2 ) * log( 2 ) +
            sum_lgamma_bnu - ( ( b + p - 1 ) / (double)2 ) * log( determinant( D, p ) );
    }
    else if(size_graph == 0)
    {
        double sum_log_diag_D = 0;
        for (int i=0; i<p; i++) sum_log_diag_D += log(D[i*(p+1)]);
        result = ( (double)p * b / 2 ) * log( 2 ) + p * lgamma( (double)b / 2 ) - ( (double)b / 2 ) * sum_log_diag_D;
    }
    else
    {
        vector<double> f_T(iter);
        log_exp_mc(G, &nu[0], b, &H[0], check_H, iter, p, &f_T[0]);
        double log_Ef_T = 0;
        for (int i=0; i<iter; i++) log_Ef_T += exp(-f_T[i]/2)/iter;
        log_Ef_T = log(log_Ef_T);
        double sum_lgamma_bnu = 0;
        for (int i=0; i<p; i++) sum_lgamma_bnu += lgamma((b + nu[i])/(double)2);
        double sum_bnu_log_diag_Ti = 0;
        vector<double> colsums(p,0);
        for (int j=0; j<p; j++) for (int i=0; i<=j; i++) colsums[j] += G[i*p + j];
        for (int i=0; i<p; i++) sum_bnu_log_diag_Ti += (b + nu[i] + colsums[i])*log(Ti[i*(p+1)]);
        double result = ( (double)size_graph / 2 ) * log( M_PI ) + ( (double)p * b / 2 + size_graph ) * log( 2 ) +
    	    sum_lgamma_bnu + sum_bnu_log_diag_Ti;
        result += log_Ef_T;
    }
    return(result);
}

int update_b(int b_old, int G[], double K[], double D[], int p, int L, int ns[])
{
    double lambda = 1 ; //What is lambda
    double q_ratio; // to be used in prior
    int b_new;
    if (b_old == 3)
    {
        b_new = 4;
        q_ratio = 0.5;
    }
    else
    {
        if ((double)rand()/RAND_MAX > 0.5) b_new = b_old + 1;
        else b_new = b_old - 1;
        if(b_old == 4 && b_new == 3) q_ratio = 2; // the case if b goes from 4 to 3.
        else q_ratio = 1;
    }
    //prior
    double logr = log(lambda)*(b_new - b_old) + log((double)factorial(b_old-3)/factorial(b_new-3)) + log(q_ratio);
    for (int l=0; l<L; l++)
    {
        logr += ns[l]*(gnorm(&G[l*p*p], b_old, D, p, 100)-gnorm(&G[l*p*p], b_new, D, p, 100)+log(determinant(&K[l*p*p], p))*(b_new - b_old));
    }
    if (logr > log((double)rand()/RAND_MAX)) return(b_new);
    else return(b_old);
}

void update_D(double D[], int G[], double K[], int b, int p, int L, int ns[])
{
    double alpha = 2.0, beta = 2.0;
    double c = 5;
    vector<double> D_new(p*p);
    memcpy(&D_new[0], D, p*p*sizeof(double));
    default_random_engine generator;
    vector<double> gnorm_old(L);
    vector<double> gnorm_new(L);
    //Calculate the first gnorm.
    for (int l=0; l<L; l++) gnorm_old[l] = gnorm(&G[l*p*p], b, D, p, 100);
    //Begin update.
    for (int j=0; j<p; j++)
    {
        double d_old = D[j*(p+1)];
        gamma_distribution<double> rgamma( c*d_old, 1/c );
        double d_new = rgamma(generator); D_new[j*(p+1)] = d_new;
        
        // Get q ratio first
        double logr = c*log(c)*(d_new - d_old) - lgamma(c*d_new) + lgamma(c*d_old) + (c*d_new-1)*log(d_old) - (c*d_old-1)*log(d_new) +c*(d_new-d_old)
        // Then prior.
                        + (alpha-1)*log(d_new / d_old) - beta*(d_new - d_old);
        for (int l=0; l<L; l++)
        {
            
            gnorm_new[l] = gnorm(&G[l*p*p], b, &D_new[0], p, 100);
            logr += ns[l]*(gnorm_old[l] - gnorm_new[l]) +log(K[l*p*p + j*(p+1)]*(d_new-d_old)) ;
        }
        if (logr > log((double)rand()/RAND_MAX))
        {
            D[j*(p+1)] = d_new;
            memcpy(&gnorm_old[0], &gnorm_new[0], L*sizeof(double));
        }
        else
        {
            D_new[j*(p+1)] = d_old; //if the last proposal is not accepted, we reload the original value from D.
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmatric matrix A (p x p) and 
// retrieves A12(1x(p-1)) and A22((p-1)x(p-1))
// Like A12=A[j, -j], and A22=A[-j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices1( double A[], double A12[], double A22[], int psub, int pdim )
{
    int i, ixpdim, ixp1, p1 = pdim - 1, subxp = psub * pdim, mpsub = pdim - psub - 1;
    int size_psub  = sizeof( double ) * psub;
    int size_mpsub = sizeof( double ) * mpsub;

    memcpy( A12,        A + subxp,            size_psub );    
    memcpy( A12 + psub, A + subxp + psub + 1, size_mpsub );    

    for( i = 0; i < psub; i++ )
    {    
        ixpdim = i * pdim;
        ixp1   = i * p1;
        
        memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
        memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
    }

    for( i = psub + 1; i < pdim; i++ )
    {
        ixpdim = i * pdim;
        ixp1   = ( i - 1 ) * p1;
        
        memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
        memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmetric matrix A (p x p) and 
// retrieves upper part of sub_matrix B (p_sub x p_sub), dictated by vector sub
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrix_upper( double A[], double sub_A[], int sub[], int psub, int pdim )
{
    int i, j, ixp, subixp;
            
    for( i = 0; i < psub; i++ )
    {
        ixp    = i * psub;
        subixp = sub[ i ] * pdim;
            
        for( j = 0; j <= i; j++ )
            sub_A[ ixp + j ] = A[ subixp + sub[ j ] ]; 
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_row_mins( double A[], double sub_A[], int subj, int pdim )
{
    int subxp = subj * pdim;

    memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );        
    memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(p-2 x 2) which is sub cols of matrix A, minus two elements
// Likes A[-(i,j), (i,j)] in R 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_cols_mins( double A[], double sub_A[], int subi, int subj, int pdim )
{    
    int p2 = pdim - 2, subixp = subi * pdim, subjxp = subj * pdim;

    memcpy( sub_A           , A + subixp           , sizeof( double ) * subi );        
    memcpy( sub_A + subi    , A + subixp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );    
    memcpy( sub_A + subj - 1, A + subixp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );    

    memcpy( sub_A + p2           , A + subjxp           , sizeof( double ) * subi );        
    memcpy( sub_A + p2 + subi    , A + subjxp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );    
    memcpy( sub_A + p2 + subj - 1, A + subjxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves A11_inv ( 2 x 2 ), A21 ( ( p - 2 ) x 2 ), and A22 ( ( p - 2 ) x ( p - 2 ) )
// Like A11_inv=inv ( A[ e, e ] ), A21 = A[ -e, e ], and A22 = A[ -e, -e ] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int sub0, int sub1, int pdim )
{
    int i, ixp, ixp2, p2 = pdim - 2;
    int sub0xp = sub0 * pdim, sub1xp = sub1 * pdim, sub0_plus = sub0 + 1, sub1_plus = sub1 + 1;
    
    double a11 = A[ sub0 * pdim + sub0 ];
    double a12 = A[ sub0 * pdim + sub1 ];
    double a22 = A[ sub1 * pdim + sub1 ];

    double det_A11 = a11 * a22 - a12 * a12;
    A11_inv[ 0 ]   = a22 / det_A11;
    A11_inv[ 1 ]   = - a12 / det_A11;
    A11_inv[ 2 ]   = A11_inv[ 1 ];
    A11_inv[ 3 ]   = a11 / det_A11;
    
    int size_sub0      = sizeof( double ) * sub0;
    int size_sub1_sub0 = sizeof( double ) * ( sub1 - sub0_plus );
    int size_pdim_sub0 = sizeof( double ) * ( pdim - sub1_plus );
    
    memcpy( A21           , A + sub0xp            , size_sub0 );        
    memcpy( A21 + sub0    , A + sub0xp + sub0_plus, size_sub1_sub0 );    
    memcpy( A21 + sub1 - 1, A + sub0xp + sub1_plus, size_pdim_sub0 );    

    memcpy( A21 + p2           , A + sub1xp            , size_sub0 );        
    memcpy( A21 + p2 + sub0    , A + sub1xp + sub0_plus, size_sub1_sub0 );    
    memcpy( A21 + p2 + sub1 - 1, A + sub1xp + sub1_plus, size_pdim_sub0 );    
 
    for( i = 0; i < sub0; i++ )
    {    
        ixp  = i * pdim;
        ixp2 = i * p2;

        memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
        memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
        memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );    
    }
 
    for( i = sub0_plus; i < sub1; i++ )
    {
        ixp  = i * pdim;
        ixp2 = ( i - 1 ) * p2;

        memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
        memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
        memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );    
    }
    
    for( i = sub1_plus; i < pdim; i++ )
    {
        ixp  = i * pdim;
        ixp2 = ( i - 2 ) * p2;
                
        memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
        memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
        memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );        
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// computing birth/death rate or alpha for element (i,j)
// it is for double Metropolis-Hasting algorihtms
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_H_ij( double K[], double sigma[], double *log_Hij, int selected_edge_i, int selected_edge_j,
               double Kj12[], double Kj12xK22_inv[], double K12[], double K12xK22_inv[], double K121[], 
               double sigmaj12[], double sigmaj22[], double sigma12[], double sigma22[], double sigma11_inv[], double sigma21xsigma11_inv[],
               int dim, int p1, int p2, int jj,
               double Dsijj, double Dsij, double Dsjj )
{
    int one = 1, two = 2;
    double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
    char transT = 'T', transN = 'N', sideL = 'L';                                                                    
    
    //double sigmaj11 = sigma[*jj];        // sigma[j, j]  
    sub_matrices1( sigma, sigmaj12, sigmaj22, selected_edge_j, dim );

    // sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
    // Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
    double sigmajj_inv = - 1.0 / sigma[ selected_edge_j * ( dim + 1 ) ];
    //F77_NAME(dsyr)( &sideL, p1, &sigmajj_inv, sigmaj12, &one, sigmaj22, p1 FCONE );
    dsyr_( &sideL, &p1, &sigmajj_inv, sigmaj12, &one, sigmaj22, &p1 );

    // For (i,j) = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |    
    sub_row_mins( K, Kj12, selected_edge_j, dim );  // K12 = K[j, -j]  
    Kj12[ selected_edge_i ] = 0.0;                       // K12[1,i] = 0

    // Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
    //F77_NAME(dsymv)( &sideL, p1, &alpha, &sigmaj22[0], p1, Kj12, &one, &beta, Kj12xK22_inv, &one FCONE );
    dsymv_( &sideL, &p1, &alpha, &sigmaj22[0], &p1, Kj12, &one, &beta, Kj12xK22_inv, &one );
    
    // K022 = Kj12xK22_inv %*% t(Kj12)
    //double K022 = F77_NAME(ddot)( p1, Kj12xK22_inv, &one, Kj12, &one );
    double K022 = ddot_( &p1, Kj12xK22_inv, &one, Kj12, &one );

    // For (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
    sub_cols_mins( K, K12, selected_edge_i, selected_edge_j, dim );   // K21 = K[-e, e] 
    
    sub_matrices_inv( sigma, sigma11_inv, sigma12, sigma22, selected_edge_i, selected_edge_j, dim );

    // sigma21xsigma11_inv = sigma21 %*% sigma11_inv
    //F77_NAME(dgemm)( &transN, &transN, p2, &two, &two, &alpha, sigma12, p2, sigma11_inv, &two, &beta, sigma21xsigma11_inv, p2 FCONE FCONE );
    dgemm_( &transN, &transN, &p2, &two, &two, &alpha, sigma12, &p2, sigma11_inv, &two, &beta, sigma21xsigma11_inv, &p2 );
    
    // sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
    //F77_NAME(dgemm)( &transN, &transT, p2, p2, &two, &alpha1, sigma21xsigma11_inv, p2, sigma12, p2, &beta1, sigma22, p2 FCONE FCONE );
    dgemm_( &transN, &transT, &p2, &p2, &two, &alpha1, sigma21xsigma11_inv, &p2, sigma12, &p2, &beta1, sigma22, &p2 );

    // K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
    //F77_NAME(dgemm)( &transT, &transN, &two, p2, p2, &alpha, K12, p2, sigma22, p2, &beta, K12xK22_inv, &two FCONE FCONE );  
    dgemm_( &transT, &transN, &two, &p2, &p2, &alpha, K12, &p2, sigma22, &p2, &beta, K12xK22_inv, &two );
    
    // K121 = K12xK22_inv %*% K21                                                    
    //F77_NAME(dgemm)( &transN, &transN, &two, &two, p2, &alpha, K12xK22_inv, &two, K12, p2, &beta, K121, &two FCONE FCONE );
    dgemm_( &transN, &transN, &two, &two, &p2, &alpha, K12xK22_inv, &two, K12, &p2, &beta, K121, &two );
    
    // Finished (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

    double a11      = K[selected_edge_i *dim + selected_edge_i] - K121[0];    
    double sum_diag = Dsjj * ( K022 - K121[3] ) - Dsij * ( K121[1] + K121[2] );

    // Dsijj = Dsii - Dsij * Dsij / Dsjj;
    //*log_Hij = ( log( static_cast<double>(*Dsjj) ) - log( static_cast<double>(a11) ) + *Dsijj * a11 - sum_diag ) / 2;
    *log_Hij = 0.5 * ( log( Dsjj / a11 ) + Dsijj * a11 - sum_diag );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// rgwish ONLY for inside of MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rgwish_sigma( int G[], int size_node[], double Ts[], double K[], double sigma[], int b_star, 
                    int dim, double threshold,
                    double sigma_start[], double inv_C[], double beta_star[], double sigma_i[], 
                    double sigma_start_N_i[], double sigma_N_i[], int N_i[] )
{
    int i, i1, j, ij, ip, l, size_node_i, info, one = 1, pxp = dim * dim, dim1 = dim + 1, bKdim = b_star + dim - 1;    
    
    double alpha = 1.0, beta  = 0.0;    
    
    char transT  = 'T', transN = 'N', side = 'R', upper = 'U';                                                                    
    
    // - - STEP 1: sampling from wishart distributions  - - - - - - - - - - - - - - - - - - - - - -|
    // - -  Sample values in Psi matrix
    //GetRNGstate();
    default_random_engine generator;
    normal_distribution<double> distribution_2(0.0,1.0);
    #pragma omp parallel for
    for( i = 0; i < dim; i++ )
    {
        //sigma_start[ i * dim1 ] = sqrt( Rf_rgamma( ( bKdim - i ) * 0.5, 2.0 ) ); // i * dim1 = i * dim + i
        gamma_distribution<double> distribution(( bKdim - i ) * 0.5, 2.0 ); //Rf_rgamma uses scale while C++ uses beta.
        sigma_start[ i * dim1 ] = sqrt( distribution(generator) ); // i * dim1 = i * dim + i
        //sigma_start[i * dim1] = sqrt( rchisq( bKdim - i ) ); // i * dim1 = i * dim + i
    }
    #pragma omp parallel for
    for( j = 1; j < dim; j++ )
        for( int i = 0; i < j; i++ )
        {
            
            sigma_start[ j * dim + i ] = distribution_2(generator);
            sigma_start[ i * dim + j ] = 0.0;
        }
    //PutRNGstate();
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
    
    // C = psi %*% Ts   I used psi = psi %*% Ts   Now is  sigma_start = sigma_start %*% Ts
    //F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &sigma_start[0], &dim FCONE FCONE FCONE FCONE );
    dtrmm_(&side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &sigma_start[0], &dim);

    side = 'L';
    // creating an identity matrix
    #pragma omp parallel for
    for( i = 0; i < dim; i++ )
        for( int j = 0; j < dim; j++ )
            inv_C[ j * dim + i ] = ( i == j );    
    
    // op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
    //F77_NAME(dtrsm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, &sigma_start[0], &dim, &inv_C[0], &dim FCONE FCONE FCONE FCONE );
    dtrsm_(&side, &upper, &transN, &transN, &dim, &dim, &alpha, &sigma_start[0], &dim, &inv_C[0], &dim);
 
    // sigma_start = inv_C %*% t( inv_C )                                                                                  
    //F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &dim, &alpha, &inv_C[0], &dim, &inv_C[0], &dim, &beta, &sigma_start[0], &dim FCONE FCONE );
    dgemm_(&transN, &transT, &dim, &dim, &dim, &alpha, &inv_C[0], &dim, &inv_C[0], &dim, &beta, &sigma_start[0], &dim);
    
    memcpy( sigma, &sigma_start[0], sizeof( double ) * pxp ); 
    
//    double temp, max_diff = 1.0, threshold_c = *threshold;
    double mean_diff = 1.0;
    int counter = 0;
    while( ( mean_diff > threshold ) and ( counter < 5000 ) )
    {
        counter++;
        mean_diff = 0.0;
        
        for( i = 0; i < dim; i++ )
        {
            ip = i * dim;

            size_node_i = size_node[ i ];
            if( size_node_i > 0 )
            {
                l = 0;
                for( j = 0; j < dim; j++ )
                {
                    ij = ip + j;
                    if( G[ ij ] )
                    {
                        sigma_start_N_i[ l ] = sigma_start[ ij ]; 
                        N_i[ l++ ]           = j;
                    }
                    else
                        beta_star[ j ] = 0.0; 
                }
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
                
                sub_matrix_upper( sigma, &sigma_N_i[0], &N_i[0], size_node_i, dim );
                    
                // A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
                //F77_NAME(dposv)( &upper, &size_node_i, &one, &sigma_N_i[0], &size_node_i, &sigma_start_N_i[0], &size_node_i, &info FCONE );
                dposv_(&upper, &size_node_i, &one, &sigma_N_i[0], &size_node_i, &sigma_start_N_i[0], &size_node_i, &info);

                for( j = 0; j < size_node_i; j++ ) beta_star[ N_i[ j ] ] = sigma_start_N_i[ j ];
    
                // sigma_i = sigma %*% beta_star
                //F77_NAME(dsymv)( &side, &dim, &alpha, sigma, &dim, &beta_star[0], &one, &beta, &sigma_i[0], &one FCONE );
                dsymv_(&side, &dim, &alpha, sigma, &dim, &beta_star[0], &one, &beta, &sigma_i[0], &one);
                
                memcpy( sigma + ip, sigma_i, sizeof( double ) * i );    
                
                for( j = 0; j < i; j++ )
                {
                    ij         = j * dim + i;
                    mean_diff += fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//                    temp      = fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//                    max_diff  = ( temp > max_diff ) ? temp : max_diff;                     

                    sigma[ ij ] = sigma_i[ j ];
                }
                
                i1 = i + 1;
                memcpy( sigma + ip + i1, sigma_i + i1, sizeof( double ) * ( dim - i1 ) );    

                for( j = i1; j < dim; j++ )
                {
                    ij         = j * dim + i;
                    mean_diff += fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//                    temp      = fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//                    max_diff  = ( temp > max_diff ) ? temp : max_diff;                     

                    sigma[ ij ] = sigma_i[ j ];
                }
            }else{                    
                for( j = 0; j < i; j++ )
                {
                    ij         = j * dim + i;
                    mean_diff += fabs( static_cast<double>( sigma[ ij ] ) );
//                    temp     = fabs( static_cast<double>( sigma[ ij ] ) );
//                    max_diff = ( temp > max_diff ) ? temp : max_diff;                     

                    sigma[ ij ]     = 0.0;
                    sigma[ ip + j ] = 0.0;
                }
                
                for( j = i + 1; j < dim; j++ )
                {
                    ij         = j * dim + i;
                    mean_diff += fabs( static_cast<double>( sigma[ ij ] ) );
//                    temp     = fabs( static_cast<double>( sigma[ ij ] ) );
//                    max_diff = ( temp > max_diff ) ? temp : max_diff;                     

                    sigma[ ij     ] = 0.0;
                    sigma[ ip + j ] = 0.0;                
                }
            } 
        }
        
        mean_diff /= pxp;
    }
    
    memcpy( &sigma_start[0], sigma, sizeof( double ) * pxp );         
    
    inverse( &sigma_start[0], K, dim );
    // creating an identity matrix
    //#pragma omp parallel for
    //for( i = 0; i < dim; i++ )
    //    for( int j = 0; j < dim; j++ )
    //        K[ j * dim + i ] = ( i == j );
    
    // LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
    //F77_NAME(dposv)( &upper, &dim, &dim, &sigma_start[0], &dim, K, &dim, &info );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel Computation for birth-death rates for double BD-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rates_bdmcmc_dmh_parallel( double rates[], double log_ratio_g_prior[], int G[], int index_row[], int index_col[], int *sub_qp, double Ds[], double D[],
                            double sigma[], double K[], double sigma_dmh[], 
                            double K_dmh[], int L, int dim )
{
    int p1 = dim - 1, p2 = dim - 2, p2x2 = ( dim - 2 ) * 2;

    #pragma omp parallel
    {
        int pxp = dim*dim;
        int index_rate_j, i, j, ij, jj;
        double Dsjj, Dsij, Dsijj, Dij, Dijj, Djj, log_rate;

        double *K121                = new double[ 4 ];  
        double *Kj12                = new double[ p1 ];  
        double *sigmaj12            = new double[ p1 ];  
        double *sigmaj22            = new double[ p1 * p1 ];  
        double *Kj12xK22_inv        = new double[ p1 ];  
        double *K21                 = new double[ p2x2 ];  
        double *sigma12             = new double[ p2x2 ];  
        double *sigma22             = new double[ p2 * p2 ];  
        double *sigma11_inv         = new double[ 4 ];  
        double *sigma21xsigma11_inv = new double[ p2x2 ];  
        double *K12xK22_inv         = new double[ p2x2 ];
        
        double *K12                 = new double[ p2x2 ];
        
        for (int l=0; l<L; l++)
        {
            #pragma omp for
            for( j = 1; j < dim; j++ )
            {            
                index_rate_j = l*dim*(dim - 1)/2 + ( j * ( j - 1 ) ) / 2;

                jj   = j * dim + j;
                Dsjj = Ds[ l*dim*dim + jj ];
                Djj  = D[ jj ];

                for( i = 0; i < j; i++ )
                {
                    ij    = j * dim + i;
                    int lij = ij + l*dim*dim;
                    Dsij  = Ds[ lij ];
                    Dsijj = - Dsij * Dsij / Dsjj;
                    Dij   = D[ ij ];
                    Dijj  = - Dij * Dij / Djj;

                    double logH_ij, logI_p;
                    
                    log_H_ij( &K[l*pxp], &sigma[l*pxp], &logH_ij, i, j,
                           &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
                           &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
                           dim, p1, p2, jj,
                           Dsijj, Dsij, Dsjj );

                    log_H_ij( &K_dmh[l*pxp], &sigma_dmh[l*pxp], &logI_p, i, j,
                           &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
                           &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
                           dim, p1, p2, jj,
                           Dijj, Dij, Djj );
                    
                    //log_rate = ( G[ ij ] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );                
                    log_rate = ( G[ lij ] ) ? ( logH_ij - logI_p ) - log_ratio_g_prior[ lij ] : ( logI_p - logH_ij ) + log_ratio_g_prior[ lij ];                
                    rates[ index_rate_j + i ] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
                }
            }
        }
        delete[] K121;  
        delete[] Kj12;  
        delete[] sigmaj12;  
        delete[] sigmaj22;  
        delete[] Kj12xK22_inv;          
        delete[] K21;  
        delete[] sigma12;  
        delete[] sigma22;  
        delete[] sigma11_inv;  
        delete[] sigma21xsigma11_inv;  
        delete[] K12xK22_inv;  
        delete[] K12;  
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To select an edge for BDMCMC algorithm  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_edge( double rates[], int *index_selected_edge, double *sum_rates, int qp_star )
{
    // rates = sum_sort_rates
    vector<double>cumulative_rates( qp_star, 0.0 );
    cumulative_rates[ 0 ] = rates[ 0 ];
    for( int i = 1; i < qp_star; i++ )
        cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
    
    *sum_rates = cumulative_rates[ qp_star - 1 ];
    
    // GetRNGstate();
    //double random_value = *sum_rates * unif_rand(); // Rf_runif( 0.0, *sum_rates );
    double random_value = *sum_rates * (double)rand()/RAND_MAX; // Rf_runif( 0.0, *sum_rates );
    // PutRNGstate();

    //int counter = 0;
    //while( random_value > cumulative_rates[ counter ] )    ++counter;
    //*index_selected_edge = counter;
     
    // To start, find the subscript of the middle position.
    int lower_bound = 0;
    int upper_bound = qp_star - 1;
    int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

    while( upper_bound - lower_bound > 1 )
    {
         //if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
        ( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
        
        position = ( lower_bound + upper_bound ) / 2;
    }
    
    *index_selected_edge = ( cumulative_rates[ position ] < random_value ) ? ++position : position;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian Graphical models  
// Based on Double Metropolis-Hastings
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_DMH_bdmcmc_ma( int iteration, int burn_in, int G[], double g_prior[],
                        double K[], int dim, double threshold, double K_hat[], double p_links[],
                        int b1, int ns[], double S_dot[], double D[], int L, int print_c,
                        double beta0[], double beta1[], double X[], double Y_bar[])
{
    int index_selected_edge, selected_edge_l, selected_edge_i, selected_edge_j, selected_edge_ij;
    int ip, i, j, ij, pxp = dim * dim, Lxpxp = L * pxp, one = 1;
    int qp = L * dim * ( dim - 1 ) / 2;
    double sum_weights = 0.0, weight_C, sum_rates;
        
    vector<double> sigma( Lxpxp ); 
    for(int l=0; l<L; l++)
        inverse( &K[l*pxp], &sigma[l*pxp], dim );            
    
    vector<double> p_links_Cpp( Lxpxp, 0.0 ); 
    vector<double> K_hat_Cpp( Lxpxp, 0.0 ); 
    // - -  for rgwish_sigma 
    vector<double> sigma_start( Lxpxp ); 
    vector<double> inv_C( Lxpxp ); 
    vector<double> beta_star( L*dim ); 
    vector<double> sigma_i( L*dim ); 
    vector<double> sigma_start_N_i( L*dim );   // For dynamic memory used
    vector<double> sigma_N_i( Lxpxp );         // For dynamic memory used
    vector<int> N_i( L*dim );                  // For dynamic memory used
    // - - - - - - - - - - - - - - - - - - - - 
    vector<double> sigma_dmh( Lxpxp );          // for double Metropolis-Hastings
    vector<double> K_dmh( Lxpxp );              // for double Metropolis-Hastings
    
    vector<double> S(Lxpxp);
    vector<double> alpha0(dim); //The intercept for mean, defaults to 0
    vector<double> alpha1(dim); //The slope for mean, defaults to 0
    vector<double> mu(L*dim); // Stores mu to save time
    
    vector<double> Ds(Lxpxp); //Newly introduced to update with every D.
    vector<double> Ts(Lxpxp); //Newly introduced to update with every D.
    vector<double> Ti(pxp); //Newly introduced to update with every D.
    
    // Counting size of notes
    vector<int> size_node( L*dim, 0 );
    for( i = 0; i < L*dim; i++ )
    {
        ip = i * dim;
        for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
    }
    
    // For finding the index of rates
    vector<int> index_L( qp ); // newly added to record group
    vector<int> index_row( qp ); // already multiplied by L above.
    vector<int> index_col( qp );
    int counter = 0 ;
    for( int l=0; l<L; l++)
        for( j = 1; j < dim; j++ )
            for( i = 0; i < j; i++ )
            {
                ij = g_prior[ (l*dim + j) * dim + i ];
                if( ( ij != 0.0 ) or ( ij != 1.0 ) )
                {
                    index_L[ counter ] = l;
                    index_row[ counter ] = i;
                    index_col[ counter ] = j;
                    counter++;
                }
            }
    int sub_qp = counter;
    cout << "sub_qp is " << sub_qp << endl;
    vector<double> rates( sub_qp );
    vector<double> log_ratio_g_prior( L * pxp );
    vector<int> b_star(L);

// - - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
    //GetRNGstate();
    cout << "Begin MCMC.\n";
    int print_conter = 0;
    for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
    {
        if( ( i_mcmc + 1 ) % print_c == 0 ){
            ++print_conter;
            //( print_conter != 20 ) ? Rprintf( "%i%%->", print_conter * 5 ) : Rprintf( " done" );
            ( print_conter != 20 ) ? cout<< print_conter * 5 << "%->" : cout << " done\n" ;
        }
        
        //Calculate mu
        // not necessary, update directly in when alpha is updated.
        /*
        for (int l=0; l<L; l++)
            for (int i; i<p; i++)
                mu[l*p+i] = alpha0[i] + alpha1[i]*X[l];
        */ 
        
        // Calculate S
        for (int l=0; l<L; l++)
        {
            for (int i=0; i<dim; i++)
            {
                S[l*pxp + i*(dim+1)] = S_dot[l*pxp + i*(dim+1)] + mu[l*dim+i]*ns[l]*(mu[l*dim+i]-2*Y_bar[l*dim+i]) ;
                for (int j=0; j<i; j++)
                {
                    S[l*pxp + i*dim + j] = S_dot[l*pxp + i*dim + j]
                        + ns[l]*(mu[l*dim+i]*(mu[l*dim+j]-Y_bar[l*dim+j])-mu[l*dim+j]*Y_bar[l*dim+i]);
                    S[l*pxp + j*dim + i] = S[l*pxp + i*dim + j];
                }
            }
        }
        
        // Calculate b_star
        for (int l=0; l<L; l++) b_star[l] = ns[l] + b1;
        //Update Ds after D is updated.
        for (int l=0; l<L; l++)
            for (int i=0; i<dim; i++)
                for (int j=0; j<dim; j++)
                {
                    int pos = i*dim + j;
                    Ds[l*pxp + pos] = D[pos] + S[l*pxp + pos];
                }
        // Subsequently update Ts.
        char uplo = 'U';
        int info;
        for (int l=0; l<L; l++)
        {
            inverse(&Ds[l*pxp], &Ts[l*pxp], dim); //Cycle Ts to save memory.
            dpotrf_( &uplo, &dim, &Ts[l*pxp], &dim, &info );
        }
        
        //Update Ti
        inverse(&D[0], &Ti[0], dim); //Cycle Ti to save memory.
        dpotrf_( &uplo, &dim, &Ti[0], &dim, &info );
        
        for (int l=0; l<L; l++)
        {
            counter = 0;
            for( j = 1; j < dim; j++ )
                for( i = 0; i < j; i++ )
                {
                    ij = (l*dim + j) * dim + i;
                    //log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
                    double odds = exp(beta0[counter] + beta1[counter]*X[l]);
                    log_ratio_g_prior[ ij ] = log( odds/ (1+odds) );
                    counter++;
                }
        }
        
// - - - STEP 1: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - -|        
        //cout << "Begin rgwish_sigma." << endl;
        // sampling from K and sigma for double Metropolis-Hastings
        for (int l=0; l<L; l++)
        {
            rgwish_sigma( &G[l*pxp], &size_node[l*dim], &Ti[0], &K_dmh[l*pxp], &sigma_dmh[l*pxp], b1, dim, threshold, &sigma_start[l*pxp],
                          &inv_C[l*pxp], &beta_star[l*dim], &sigma_i[l*dim], &sigma_start_N_i[l*dim], &sigma_N_i[l*pxp], &N_i[l*dim] );
        }
        
        //cout << "Calculate rates." << endl;
        rates_bdmcmc_dmh_parallel( &rates[0], &log_ratio_g_prior[0], G, &index_row[0], &index_col[0],
                                   &sub_qp, &Ds[0], D, &sigma[0], K, &sigma_dmh[0], &K_dmh[0], L, dim );
        
        // Selecting an edge based on birth and death rates
        //cout << "Select edge." << endl;
        select_edge( &rates[0], &index_selected_edge, &sum_rates, sub_qp );
        selected_edge_l = index_L[ index_selected_edge ];
        selected_edge_i = index_row[ index_selected_edge ];
        selected_edge_j = index_col[ index_selected_edge ];

        //Update b and D.
        b1 = update_b(b1, G, K, D, dim, L, &ns[0]); // Note that this b1 will not be sent back the main.
        update_beta(beta0, beta1, G, X, dim, L, ns);
        update_D(D, G, K, b1, dim, L, &ns[0]);
        
        //Update alpha.
        update_alpha( &alpha0[0], &alpha1[0], &K[0], &X[0], &Y_bar[0], &mu[0], &ns[0], dim, L);

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|    
        if( i_mcmc >= burn_in )
        {
            weight_C = 1.0 / sum_rates;
            
            // K_hat_Cpp[i] += K[i] / sum_rates;
            //F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
            //LAPACK daxpy
            // http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga8f99d6a644d3396aa32db472e0cfc91c.html
            daxpy_(&Lxpxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one);
            
            for( i = 0; i < Lxpxp ; i++ )
                if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
            
            sum_weights += weight_C;
        } 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |    
        
        // Updating G (graph) based on selected edge
        selected_edge_ij    = (selected_edge_l*dim + selected_edge_j) * dim + selected_edge_i;
        G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
        G[ (selected_edge_l*dim + selected_edge_i) * dim + selected_edge_j ] = G[ selected_edge_ij ];

        if( G[ selected_edge_ij ] )
        { 
            ++size_node[ selected_edge_l*dim + selected_edge_i ]; 
            ++size_node[ selected_edge_l*dim + selected_edge_j ]; 
        }else{ 
            --size_node[ selected_edge_l*dim + selected_edge_i ]; 
            --size_node[ selected_edge_l*dim + selected_edge_j ]; 
        }

// - - -- STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
        for (int l=0; l<L; l++)
            rgwish_sigma( &G[l*pxp], &size_node[l*dim], &Ts[l*pxp], &K[l*pxp], &sigma[l*pxp], b_star[l], dim, threshold, &sigma_start[l*pxp],
                          &inv_C[l*pxp], &beta_star[l*dim], &sigma_i[l*dim], &sigma_start_N_i[l*dim], &sigma_N_i[l*pxp], &N_i[l*dim] );
        
    }  
    //PutRNGstate();
    
// - - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

    for( i = 0; i < Lxpxp; i++ )
    {    
        p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
        K_hat[ i ]   = K_hat_Cpp[ i ] / sum_weights;
    }
}

void read_file_int(vector<int> *values, string file_name)
{
    ifstream inputFile(file_name);        // Input file stream object
    // Check if exists and then open the file.
    if (inputFile.good()) {
        // Push items into a vector
        int current_number;
        while (inputFile >> current_number){
            (*values).push_back(current_number);
        }
    // Close the file.
        inputFile.close();
    }else {
        cout << "Error in reading "<< file_name <<"!" << endl;
        exit(1);
    }
}

void read_file_double(vector<double> *values, string file_name)
{
    ifstream inputFile(file_name);        // Input file stream object
    // Check if exists and then open the file.
    if (inputFile.good()) {
        // Push items into a vector
        double current_number;
        while (inputFile >> current_number){
            (*values).push_back(current_number);
        }
    // Close the file.
        inputFile.close();
    }else {
        cout << "Error in reading "<< file_name <<"!" << endl;
        exit(1);
    }
}

void write_file_int(vector<int> *values, string file_name, int L, int row, int col)
{
    ofstream myfile (file_name);
    if (myfile.is_open())
    {
        for (int l=0; l<L; l++)
        {
            myfile << "Group " << l+1 << ":" << endl;
            for (int i =0; i<row; i++)
            {
                for (int j=0; j<col; j++)
                {
                    myfile << (*values)[ (l*row + i)*col + j ] << "\t";
                }
                myfile << endl;
            }
        }
        myfile.close();
    }
    else cout << "Unable to open file";
    
}

void write_file_double(vector<double> *values, string file_name, int L, int row, int col)
{
    ofstream myfile (file_name);
    if (myfile.is_open())
    {
        for (int l=0; l<L; l++)
        {
            myfile << "Group " << l+1 << ":" << endl;
            for (int i =0; i<row; i++)
            {
                for (int j=0; j<col; j++)
                {
                    myfile << (*values)[ (l*row + i)*col + j ] << "\t";
                }
                myfile << endl;
            }
        }
        myfile.close();
    }
    else cout << "Unable to open file";
    
}

int main(int argc, char **argv)
{
    // 0. Read in arguments to function
    int iter = 5000;
    int burnin = 2500;
    int p = 5;
    double threshold = .00000001;
    int b = 7;
    int L = 4; // No. of groups
    unsigned int seed = 12345;
    string g_prior_file_name = "demo_data/demo_g_prior.txt";
    string K_file_name = "demo_data/demo_K.txt";
    string ns_file_name = "demo_data/demo_n.txt";
    string S_dot_file_name = "demo_data/demo_S.txt";
    string D_file_name = "demo_data/demo_D.txt";
    string beta0_file_name = "demo_data/demo_beta0.txt";
    string beta1_file_name = "demo_data/demo_beta1.txt";
    string X_file_name = "demo_data/demo_X.txt";
    string Y_bar_file_name = "demo_data/demo_Y_bar.txt";
    string out_dir = "demo_out_change_means/";
    
    int opt;
    
    while ((opt = getopt (argc, argv, "i:b:p:L:n:B:g:K:S:D:z:o:X:Y:e:s:O:")) != -1){
        switch(opt){
            case 'i':
                iter = atoi(optarg);
                break;
            case 'b':
                burnin = atoi(optarg);
                break;
            case 'p':
                p = atoi(optarg);
                break;
            case 'L':
                L = atoi(optarg);
                break;
            case 'n':
                ns_file_name = optarg;
                break;
            case 'B':
                b = atoi(optarg);
                break;
            case 'g':
                g_prior_file_name = optarg;
                break;
            case 'K':
                K_file_name = optarg;
                break;
            case 'S':
                S_dot_file_name = optarg;
                break;
            case 'D':
                D_file_name = optarg;
                break;
            case 'z':
                beta0_file_name = optarg;
                break;
            case 'o':
                beta1_file_name = optarg;
                break;
            case 'X':
                X_file_name = optarg;
                break;
            case 'Y':
                Y_bar_file_name = optarg;
                break;
            case 'e':
                threshold = atof(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'O':
                out_dir = optarg;
                break;
            default:
            cout << "Error Usage: "<< argv[0] << "[-i (i)terations] [-b (b)urnin] [-p dimension of Y] [-L no. of groups] "
            << "[-n file of sample size of groups] [-B initial value of b] [-g file of initial g_prior] "
            << "[-K file of initial K] [-S file of S_dot] "
            << "[-D file of initial D] [-z file of initial beta_(z)ero] [-o file of initial beta_(o)ne] "
            << "[-X file of X] [-Y file of Y_bar] [-e threshold] [-s seed] [-O output directory]\n";
            exit(1);
        }
    }
    
    //Run checking on user input.
    if(burnin > iter){
        printf("Burn-in iterations must be less than total number of iterations.\n");
        return (-1);
    }
    
    int print = iter/20;
    
    // Further declare vairables
    vector<int> G(L*p*p); //L*p*p
    vector<double> g_prior; //L*p*p
    vector<double> K; //L*p*p
    vector<double> K_hat(L*p*p); //L*p*p
    vector<double> p_links(L*p*p); //L*p*p
    vector<int> ns; //L
    vector<double> S_dot; //L*p*p
    vector<double> D; //p*p
    vector<double> beta0; // p*(p-1)/2
    vector<double> beta1; // p*(p-1)/2
    vector<double> X; // L
    vector<double> Y_bar; // L*p
    
    // Read Files.
    read_file_double(&g_prior, g_prior_file_name);
    read_file_double(&K, K_file_name);
    read_file_int(&ns, ns_file_name);
    read_file_double(&S_dot, S_dot_file_name);
    read_file_double(&D, D_file_name);
    read_file_double(&beta0, beta0_file_name);
    read_file_double(&beta1, beta1_file_name);
    read_file_double(&X, X_file_name);
    read_file_double(&Y_bar, Y_bar_file_name);
    
    //for (int i=0; i<p*p; i++) cout << D[i] << endl;
    cout << "Files read." << endl;
    
    srand(seed);
    
    ggm_DMH_bdmcmc_ma(iter, burnin, &G[0], &g_prior[0], &K[0], p,
                        threshold, &K_hat[0], &p_links[0], b, &ns[0], &S_dot[0],
                        &D[0], L, print, &beta0[0], &beta1[0], &X[0], &Y_bar[0]);
    
    // Write files.
    write_file_int(&G, out_dir+"BRUG_out_G.txt",L,p,p);
    write_file_double(&p_links, out_dir+"BRUG_out_p_links.txt",L,p,p);
    write_file_double(&K, out_dir+"BRUG_out_K.txt",L,p,p);
    write_file_double(&K_hat, out_dir+"BRUG_out_K_hat.txt",L,p,p);
    
}