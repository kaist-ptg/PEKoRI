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





#ifndef LRI_OSZPROB_H
#define LRI_OSZPROB_H 1

// Definitions of internal representations for new parameters
#define LRI_ALPHAprime 6
#define LRI_ETA 7
#define LRI_BET 8
#define LRI_GAM 9
#define LRI_DEL 10



// Total number of parameters
#define LRI_TOT_NO 11


#define LRI_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
#define LRI_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#define LRI_Ne_CORE         0.468      /*   of MSW potentials                        */


// LRI potential

#define m_nucleon  1.671e-27 /* kg */
#define M_sun  1.989e30 /* kg */
#define M_earth  5.972e24 /* kg */

#define d_earthsun  7.570e17 /* /eV */ /* = 1.496e11 m */ 
#define R_earth  3.224e13 /* /eV */ /* = 6.371e6 m */ 

#define Y_sun_p  6 /* p to n ratio in the sun  */
#define Y_sun_e  6 /* e to n ratio in the sun  */
#define Y_earth_p  1 /* p to n ratio in the earth  */
#define Y_earth_e  1 /* e to n ratio in the earth  */



#include <complex.h>
#include <globes/globes.h>

// Cover functions for letting phases be free/fixed
void lri_fixAllLRI(glb_projection p);
void lri_freeAllLRI(glb_projection p);
void lri_fixEps(glb_projection p,int n);
void lri_freeEps(glb_projection p,int n);

// Prior function for extra parameters (currently equivalent to built-in priors)
double lri_prior(glb_params p,void* user_data);

// Cover function to set all LRI 
void lri_setAllLRI(glb_params p, double v);

// Modified functions originially defined in glb_probability
int lri_init_probability_engine();
int lri_free_probability_engine();
int lri_set_oscillation_parameters(glb_params p, void *user_data);
int lri_get_oscillation_parameters(glb_params p, void *user_data);
int lri_hamiltonian_cd(double E, double V, int cp_sign);
int lri_S_matrix_cd(double E, double L, double V, int cp_sign);

// Filtered probability not implemented
//int lri_filtered_probability_matrix_cd(double P[3][3], double E, double L,
//                                       double V, double sigma, int cp_sign);
int lri_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data);

// Troubleshooting
void lri_print_stored_matrices();
void lri_print_stored_N();

#endif /* LRI_OSZPROB_H */

