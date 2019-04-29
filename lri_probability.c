// PELoRI - Probability Engine for LOng-Range Interactions
// Also see: http://lotr.wikia.com/wiki/Pelori

// Modified for long-range interactions
// Sushant K. Raut
// IBS CTPU, Daejeon
// 20171204




/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <globes/globes.h>
#include "lri_probability.h"


/* Macros */
#define SQR(x)      ((x)*(x))                        /* x^2   */
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */

/* Fundamental oscillation parameters */
static double th12, th13, th23; // Mixing angles
static double delta;            // Dirac CP phase
static double mq[3];            // Squared masses

/* PELoRI parameters */
static double eps[3][3];        // LRI potential matrix
static double leta=0, lbet=0, lgam=0, ldel=0, lalphaprime=0;

// Y parameter
static double Y = LRI_Ne_MANTLE; // Ne/NN

// LRI potential
double Vnew = 0.0; // defined in the lri_init_probability_engine() function
double p0=0, p1=0, p2=0, p3=0;


/* Internal temporary variables - storing everything in mass basis */
gsl_matrix_complex *lri_U=NULL; /* The vacuum mixing matrix                           */
gsl_matrix_complex *lri_H=NULL; /* Neutrino Hamiltonian                               */
gsl_matrix_complex *lri_Q=NULL; /* Eigenvectors of Hamiltonian in mass basis          */
gsl_vector *lri_lambda=NULL;    /* Eigenvalues of Hamiltonian                         */
gsl_matrix_complex *lri_S=NULL; /* The neutrino S-matrix                              */


gsl_matrix_complex *lri_V_template=NULL;   /* Used in the construction of the matter Hamiltonian */
gsl_matrix_complex *lri_H0_template=NULL;  /* Used in the construction of the vac. Hamiltonian   */
gsl_matrix_complex *lri_S1=NULL, *lri_T0=NULL; /* Temporary matrix storage                           */

gsl_matrix_complex *lri_Vnew_template=NULL;   /* Used in the construction of the LRI matter Hamiltonian */



/***************************************************************************
 *                  I N T E R N A L   F U N C T I O N S                    *
 ***************************************************************************/



/***************************************************************************
 * Function lri_init_probability_engine                                    *
 ***************************************************************************
 * Allocates internal data structures for the probability engine.          *
 ***************************************************************************/
int lri_init_probability_engine()
{
  lri_free_probability_engine();
  
  lri_U = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_H = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_Q = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_lambda = gsl_vector_alloc(GLB_NU_FLAVOURS);
  lri_S = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    
  lri_V_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_H0_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_S1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lri_T0 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);

  lri_Vnew_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);


  return 0;
}


/***************************************************************************
 * Function lri_free_probability_engine                                    *
 ***************************************************************************
 * Destroys internal data structures of the probability engine.            *
 ***************************************************************************/
int lri_free_probability_engine()
{
  if (lri_T0!=NULL)     { gsl_matrix_complex_free(lri_T0);  lri_T0 = NULL; }
  if (lri_S1!=NULL)     { gsl_matrix_complex_free(lri_S1);  lri_S1 = NULL; }
  if (lri_H0_template!=NULL) { gsl_matrix_complex_free(lri_H0_template);  lri_H0_template = NULL; }
  if (lri_V_template!=NULL) { gsl_matrix_complex_free(lri_V_template); lri_V_template = NULL; }
  if (lri_Vnew_template!=NULL) { gsl_matrix_complex_free(lri_Vnew_template); lri_Vnew_template = NULL; }

  if (lri_S!=NULL)      { gsl_matrix_complex_free(lri_S);   lri_S = NULL; }
  if (lri_lambda!=NULL) { gsl_vector_free(lri_lambda);      lri_lambda = NULL; }
  if (lri_Q!=NULL)      { gsl_matrix_complex_free(lri_Q);   lri_Q = NULL; }
  if (lri_H!=NULL)      { gsl_matrix_complex_free(lri_H);   lri_H = NULL; }
  if (lri_U!=NULL)      { gsl_matrix_complex_free(lri_U);   lri_U = NULL; }

  return 0;
}

// Troubleshooting function
void lri_print_stored_matrices(){
  
  int n, m;
  printf("Vnew:\n");
  for (n = 0; n<3; n++){
    for (m = 0; m<3; m++){
      printf("[%-+12.3g]",eps[n][m]);
    }
    printf("\n");
  }
  printf("\n\n");
  printf("%g\n",Vnew);
  printf("%g %g %g %g\n",p0,p1,p2,p3);
}

// Another troubleshooting function
void lri_print_stored_N(gsl_matrix_complex *matrix){
  int n, m;
  printf("N =\n");
  double complex (*_mymatrix)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(matrix, 0, 0);
  for (n = 0; n<3; n++){
    printf("[%-+12.3g%-+12.3g%-+12.3g]\n",creal(_mymatrix[n][0]),creal(_mymatrix[n][1]),creal(_mymatrix[n][2]));
  }
  printf("+ I *\n");
  for (n = 0; n<3; n++){
    printf("[%-+12.3g%-+12.3g%-+12.3g]\n",cimag(_mymatrix[n][0]),cimag(_mymatrix[n][1]),cimag(_mymatrix[n][2]));
  }
  printf("\n\n");
}


/***************************************************************************
 * Function lri_set_oscillation_parameters                                 *
 ***************************************************************************
 * Sets the fundamental oscillation parameters and precomputes the mixing  *
 * matrix and part of the Hamiltonian.                                     *
 *                                                                         *
 * Also computes the template for the matter interactions                  *
 ***************************************************************************/
int lri_set_oscillation_parameters(glb_params p, void *user_data)
{
  lri_init_probability_engine();
  double complex (*_lri_U)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_U, 0, 0);
  double complex (*_V)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_V_template, 0, 0);
  double complex (*_Vnew)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_Vnew_template, 0, 0);
  int i;

  double Vsun=0, Vearth=0;
  double Nnsun = M_sun/(m_nucleon*(1+Y_sun_p)), Nnearth = M_earth/(m_nucleon*(1+Y_earth_p));

  int n, m; // loop variables

  /* Copy parameters */
  th12  = glbGetOscParams(p,GLB_THETA_12);
  th13  = glbGetOscParams(p,GLB_THETA_13);
  th23  = glbGetOscParams(p,GLB_THETA_23);
  delta = glbGetOscParams(p,GLB_DELTA_CP);
  mq[0] = abs(glbGetOscParams(p,GLB_DM_31));
  mq[1] = abs(glbGetOscParams(p,GLB_DM_31)) + glbGetOscParams(p,GLB_DM_21);
  mq[2] = abs(glbGetOscParams(p,GLB_DM_31)) + glbGetOscParams(p,GLB_DM_31);

  lalphaprime = glbGetOscParams(p,LRI_ALPHAprime);
  leta = glbGetOscParams(p,LRI_ETA);
  lbet = glbGetOscParams(p,LRI_BET);
  lgam = glbGetOscParams(p,LRI_GAM);
  ldel = glbGetOscParams(p,LRI_DEL);

  p0 = leta;
  p1 = -leta + lbet - ldel;
  p2 = -leta - lbet + lgam;
  p3 = -leta - lgam + ldel;


  Vsun = lalphaprime*Nnsun/d_earthsun * (p0 + Y_sun_p*p0 + Y_sun_e*p1);
  Vearth = lalphaprime*Nnearth/R_earth * (p0 + Y_earth_p*p0 + Y_earth_e*p1);
  Vnew = Vsun + Vearth;

  for(n = 0; n < 3; n++)
    for(m = 0; m < 3; m++)
      eps[n][m] = 0.0;
  eps[0][0] = p1 * Vnew;
  eps[1][1] = p2 * Vnew;
  eps[2][2] = p3 * Vnew;
   
  /* Compute standard vacuum mixing matrix */
  _lri_U[0][0] = cos(th12)*cos(th13);
  _lri_U[0][1] = sin(th12)*cos(th13);
  _lri_U[0][2] = sin(th13) * cexp(-I * delta);

  _lri_U[1][0] = -sin(th12)*cos(th23) - cos(th12)*sin(th23)*sin(th13) * cexp(I*delta);
  _lri_U[1][1] =  cos(th12)*cos(th23) - sin(th12)*sin(th23)*sin(th13) * cexp(I*delta);
  _lri_U[1][2] =  sin(th23)*cos(th13);

  _lri_U[2][0] =  sin(th12)*sin(th23) - cos(th12)*cos(th23)*sin(th13) * cexp(I*delta);
  _lri_U[2][1] = -cos(th12)*sin(th23) - sin(th12)*cos(th23)*sin(th13) * cexp(I*delta);
  _lri_U[2][2] =  cos(th23)*cos(th13);

  /* Calculate energy independent matrix H0 * E (mass basis) */
  gsl_matrix_complex_set_zero(lri_H0_template);
  gsl_matrix_complex_set_zero(lri_H);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    gsl_matrix_complex_set(lri_H0_template, i, i, gsl_complex_rect(0.5*mq[i], 0.0));

  /* Define matter interaction template (flavor basis) */
  _V[0][0] = 1.0;
  _V[0][1] = 0.0;
  _V[0][2] = 0.0;
  
  _V[1][0] = 0.0;
  _V[1][1] = 0.0;
  _V[1][2] = 0.0;
  
  _V[2][0] = 0.0;
  _V[2][1] = 0.0;
  _V[2][2] = 0.0;

  for(n = 0; n < 3; n++)
    for(m = 0; m < 3; m++)
      _Vnew[n][m] = eps[n][m];


  /* Compute interaction template (mass basis) */
  // Produce U^dagger.V0 in T0 (T0 = U^dagger.V0)
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_U, lri_V_template, GSL_COMPLEX_ZERO, lri_T0);
  // Produce U^dagger.V0.U in lri_V_template (lri_V_template = T0.U = U^dagger.V0.U)
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_T0, lri_U, GSL_COMPLEX_ZERO, lri_V_template);

  /* Compute new interaction template (mass basis) */
  // Produce U^dagger.Vnew in T0 (T0 = U^dagger.Vnew)
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_U, lri_Vnew_template, GSL_COMPLEX_ZERO, lri_T0);
  // Produce U^dagger.Vnew.U in lri_Vnew_template (lri_Vnew_template = T0.U = U^dagger.Vnew.U)
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_T0, lri_U, GSL_COMPLEX_ZERO, lri_Vnew_template);
  //lri_print_stored_N(lri_Vnew_template);

  return 0;
}


/***************************************************************************
 * Function lri_get_oscillation_parameters                                 *
 ***************************************************************************
 * Returns the current set of oscillation parameters.                      *
 *                                                                         *
 * Also includes the LRI parameters.                                       *
 ***************************************************************************/
int lri_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbDefineParams(p, th12, th13, th23, delta, mq[1] - mq[0], mq[2] - mq[0]);

  glbSetOscParams(p,lalphaprime,6);
  glbSetOscParams(p,leta,7);
  glbSetOscParams(p,lbet,8);
  glbSetOscParams(p,lgam,9);
  glbSetOscParams(p,ldel,10);

  return 0;
}



/***************************************************************************
 * Function lri_hamiltonian_cd                                             *
 ***************************************************************************
 * Calculates the Hamiltonian for neutrinos (cp_sign=1) or antineutrinos   *
 * (cp_sign=-1) with energy E, propagating in matter of density V          *
 * (> 0 even for antineutrinos) and stores the result in H.                *
 *                                                                         *
  ***************************************************************************/
int lri_hamiltonian_cd(double E, double V, int cp_sign)
{
  double inv_E = 1.0 / E;
  double complex (*_lri_H)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_H, 0, 0);
  double complex (*_lri_H0_template)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_H0_template, 0, 0);
  double complex (*_lri_V_template)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_V_template, 0, 0);
  double complex (*_lri_Vnew_template)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(lri_Vnew_template, 0, 0);
  int i, j;

  if (cp_sign > 0)
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _lri_H[i][j] = _lri_H0_template[i][j] * inv_E + V*_lri_V_template[i][j] + _lri_Vnew_template[i][j];
  }
  else
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _lri_H[i][j] = conj(_lri_H0_template[i][j] * inv_E - V*_lri_V_template[i][j] - _lri_Vnew_template[i][j]); 
    /* phases and V change sign for anti-nu */
  }


  return 0;
}


/***************************************************************************
 * Function lri_S_matrix_cd                                                *
 ***************************************************************************
 * Calculates the S matrix for neutrino oscillations in matter of constant *
 * density using a fast eigenvalue solver optimized to 3x3 matrices.       *
 ***************************************************************************
 * Parameters:                                                             *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *                                                                         *
 *                                                                         *
 ***************************************************************************/
int lri_S_matrix_cd(double E, double L, double V, int cp_sign)
{
  /* Introduce some abbreviations */
  double complex (*_lri_S)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(lri_S,0,0);
  double complex (*_lri_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(lri_Q,0,0);
  double complex (*_lri_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(lri_T0,0,0);
  double *_lri_lambda = gsl_vector_ptr(lri_lambda,0);
  int status;
  int i, j, k;
  
  /* Always assume matter */
  double complex (*_lri_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(lri_H,0,0);
  
  /* Calculate neutrino lri_Hamiltonian */
  if ((status=lri_hamiltonian_cd(E, V, cp_sign)) != 0)
    return status;
  
  /* Calculate eigenvalues of Hamiltonian */
  if ((status=zheevh3(_lri_H, _lri_Q, _lri_lambda)) != 0)
    return status;
  
  /* Calculate S-Matrix in matter basis ... */
  double phase;
  gsl_matrix_complex_set_zero(lri_S);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase    = -L * _lri_lambda[i];
    _lri_S[i][i] = cos(phase) + I*sin(phase);
  } 
  
  /* ... and transform it to the mass basis */
  gsl_matrix_complex_set_zero(lri_T0);
  double complex *p = &_lri_T0[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* lri_T0 = S.Q^\dagger */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      int k;
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_lri_S[i][k])*creal(_lri_Q[j][k])+cimag(_lri_S[i][k])*cimag(_lri_Q[j][k]) )
                + I * ( cimag(_lri_S[i][k])*creal(_lri_Q[j][k])-creal(_lri_S[i][k])*cimag(_lri_Q[j][k]) );
      }
      p++;
    }
  gsl_matrix_complex_set_zero(lri_S);
  p = &_lri_S[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* S = Q.T0 */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_lri_Q[i][k])*creal(_lri_T0[k][j])-cimag(_lri_Q[i][k])*cimag(_lri_T0[k][j]) )
                + I * ( cimag(_lri_Q[i][k])*creal(_lri_T0[k][j])+creal(_lri_Q[i][k])*cimag(_lri_T0[k][j]) );
      }
      p++;
    }

  return 0;
}


// Filetering probablity function, not yet implemented

/***************************************************************************
 * Function lri_filtered_probability_matrix_cd                             *
 ***************************************************************************
 * Calculates the probability matrix for neutrino oscillations in matter   *
 * of constant density, including a low pass filter to suppress aliasing   *
 * due to very fast oscillations.                                          *
 ***************************************************************************
 * Parameters:                                                             *
 *   P: Storage buffer for the probability matrix                          *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   sigma: Width of Gaussian filter                                       *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *                                                                         *
 * NSI addition uses NSI hamiltionian.                                     *
 ***************************************************************************/
/* int nsi_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V, */
/*                                        double sigma, int cp_sign) */
/* { */
/*   /\* Introduce some abbreviations *\/ */
/*   double complex (*_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Q,0,0); */
/*   double complex (*_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0,0,0); */
/*   double *_lambda = gsl_vector_ptr(lambda,0); */
/*   int status; */
/*   int i, j, k, l; */
 
/*   /\* Assume Matter *\/ */

/*   double complex (*_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(H,0,0); */
  
/*   /\* Calculate neutrino Hamiltonian *\/ */
/*   if ((status=nsi_hamiltonian_cd(E, V, cp_sign)) != 0) */
/*     return status; */
  
/*   /\* Calculate eigenvalues of Hamiltonian *\/ */
/*   if ((status=zheevh3(_H, _Q, _lambda)) != 0) */
/*     return status; */


/*   // Calculate probability matrix (see GLoBES manual for a discussion of the algorithm) */
/*   double phase, filter_factor; */
/*   gsl_matrix_complex_set_zero(T0); */
/*   for (i=0; i < GLB_NU_FLAVOURS; i++) */
/*     for (j=i+1; j < GLB_NU_FLAVOURS; j++) */
/*     { */
/*       phase         = -L * (_lambda[i] - _lambda[j]); */
/*       filter_factor = exp(-0.5 * SQR(phase*sigma) / SQR(1.0e-9 * E)); */
/*       _T0[i][j]     = filter_factor * (cos(phase) + I*sin(phase)); */
/*     } */

/*   for (k=0; k < GLB_NU_FLAVOURS; k++) */
/*     for (l=0; l < GLB_NU_FLAVOURS; l++) */
/*     { */
/*       P[k][l] = 0.0; */
/*       for (i=0; i < GLB_NU_FLAVOURS; i++) */
/*       { */
/*         for (j=i+1; j < GLB_NU_FLAVOURS; j++) */
/*           P[k][l] += 2 * creal(_Q[k][j]*conj(_Q[l][j])*conj(_Q[k][i])*_Q[l][i]*_T0[i][j]); */
/*         P[k][l] += SQR_ABS(_Q[k][i]) * SQR_ABS(_Q[l][i]); */
/*       } */
/*     } */
    
/*   return 0; */
/* } */


/***************************************************************************
 * Function nsi_probability_matrix                                         *
 ***************************************************************************
 * Calculates the neutrino oscillation probability matrix.                 *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:       Buffer for the storage of the matrix                         *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *   E:       Neutrino energy (in GeV)                                     *
 *   psteps:  Number of layers in the matter density profile               *
 *   length:  Lengths of the layers in the matter density profile in km    *
 *   density: The matter densities in g/cm^3                               *
 *   filter_sigma: Width of low-pass filter or <0 for no filter            *
 *   user_data: int*, tells if source and detector normalization should    *
 *                    be used                                              *
 *                                                                         *
 * NSI addition uses NSI probability calculation.                          *
 ***************************************************************************/
int lri_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data)
{
  int status;
  int i, j;

  /* Convert energy to eV */
  E *= 1.0e9;
  

  // Filtering not implemented
/*   if (filter_sigma > 0.0)                     /\* With low-pass filter *\/ */
/*   { */
/*     if (psteps == 1) */
/*       lri_filtered_probability_matrix_cd(P, E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, */
/*           filter_sigma, cp_sign); */
/*     else */
/*       return -1; */
/*   } */
/*   else   */                                      /* Without low-pass filter */
  //  {
  if (psteps > 1){
    gsl_matrix_complex_set_identity(lri_S1);                                 /* S1 = 1 */
    for (i=0; i < psteps; i++){
      status = lri_S_matrix_cd(E, GLB_KM_TO_EV(length[i]), density[i]*LRI_V_FACTOR*LRI_Ne_MANTLE, cp_sign);
      if (status != 0)
	return status;
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_S, lri_S1, /* T0 = S.S1 */
		     GSL_COMPLEX_ZERO, lri_T0);
      gsl_matrix_complex_memcpy(lri_S1, lri_T0);                                 /* S1 = T0 */
    } 
    gsl_matrix_complex_memcpy(lri_S, lri_S1);                                    /* S = S1 */
  }
  else{
    status = lri_S_matrix_cd(E, GLB_KM_TO_EV(length[0]), density[0]*LRI_V_FACTOR*LRI_Ne_MANTLE, cp_sign);
    if (status != 0)
      return status;
  }
  
  // S matrix in mass basis stored in S

  // Transform to flavor basis, use T0 as temporary matrix storage
  // We now need to compute the evolution matrix in the flavor basis by rotating with N
  if (cp_sign > 0){ // Neutrinos
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, lri_S, lri_U, GSL_COMPLEX_ZERO, lri_T0); // T0 = S.U^dagger
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, lri_U, lri_T0, GSL_COMPLEX_ZERO, lri_S); // S = U.T0 = U.S.U^dagger
  }
  else { // Antineutrinos
    gsl_blas_zgemm(CblasTrans,CblasConjTrans,GSL_COMPLEX_ONE,lri_S,lri_U,GSL_COMPLEX_ZERO,lri_T0); // T0 = S^T.U^dagger
    gsl_blas_zgemm(CblasTrans,CblasTrans,GSL_COMPLEX_ONE,lri_T0,lri_U,GSL_COMPLEX_ZERO,lri_S); // S = T0^T.U^T = U^*.S.U^T
  }

  double complex (*_lri_S)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(lri_S,0,0);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    for (j=0; j < GLB_NU_FLAVOURS; j++)
      P[j][i] = SQR_ABS(_lri_S[i][j]);
  //  }

  return 0;
}

/***************************************************************************
 * Function lri_priors                                                     *
 ***************************************************************************
 * Calculates the priors for the fit.                                      *
 ***************************************************************************
 * Parameters:                                                             *
 *   in:        The input parameters.                                      *
 *   user_data: Additional user data, not used                             *
 ***************************************************************************/

double lri_prior(const glb_params in, void* user_data){
  glb_params central_values = glbAllocParams();
  glb_params input_errors   = glbAllocParams();
  glb_projection p          = glbAllocProjection();
  glbGetCentralValues(central_values);
  glbGetInputErrors(input_errors);
  glbGetProjection(p);
  int i,N;
  double pv = 0.0;
  double fitvalue,centralvalue,inputerror;

  N = glbGetNumOfOscParams();

  /* oscillation parameter priors */
  for(i=0;i < N; i++){
    if(glbGetProjectionFlag(p,i)==GLB_FREE){
      fitvalue = glbGetOscParams(in,i);
      centralvalue = glbGetOscParams(central_values,i);
      inputerror = glbGetOscParams(input_errors,i);
      if(inputerror > 1e-12){
	double x = (centralvalue-fitvalue)/inputerror;
	pv += x*x;
      }
    }
  }

  /* matter parameter priors */
  for(i=0;i<glb_num_of_exps;i++){
    if(glbGetDensityProjectionFlag(p,i)==GLB_FREE){
      fitvalue = glbGetDensityParams(in,i);
      centralvalue = 1.0;
      inputerror = glbGetDensityParams(input_errors,i);
      if(inputerror > 1e-12){
	double x = (centralvalue-fitvalue)/inputerror;
	pv += x*x;
      }
    }
  }
  return pv;
}



// Cover functions for letting phases be free/fixed

void lri_fixAllLRI(glb_projection p){
  int k; // Loop variable
  for(k=LRI_ALPHAprime; k < LRI_TOT_NO; k++)
    glbSetProjectionFlag(p,GLB_FIXED,k);
}
void lri_freeAllLRI(glb_projection p){
  int k; // Loop variable
  for(k=LRI_ALPHAprime; k < LRI_TOT_NO; k++)
    glbSetProjectionFlag(p,GLB_FREE,k);
}
void lri_fixLRI(glb_projection p,int n){
  glbSetProjectionFlag(p,GLB_FIXED,n);
}
void lri_freeLRI(glb_projection p,int n){
  glbSetProjectionFlag(p,GLB_FREE,n);
}


// Functions for setting multiple parameters

void lri_setAllLRI(glb_params p, double v){
  int k;
  for(k = LRI_ALPHAprime; k < LRI_TOT_NO; k++)
    glbSetOscParams(p,v,k);
}
