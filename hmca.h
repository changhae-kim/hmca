#ifndef HMCA_H
#define HMCA_H

#include <math.h>
#include <stdlib.h>
#include <string.h>


typedef double (*hmca_nn) (const int *indices, int n_species_0);

typedef double (*mlmc_closure) (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);


// Shared Functions

inline int hmca_sym_id (int i, int j, int n)
{
	return (i < j) ? (n*i+j-i*(i+1)/2) : (n*j+i-j*(j+1)/2);
}

void hmca_lognorm_1d (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound);

void hmca_logexp_1d (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound);


// Mean-Field Approximation

inline double hmca_mf_nn_1x1 (const int *indices, int n_species_0)
{
	return 4.0;
}

inline double hmca_mf_nn_2x1 (const int *indices, int n_species_0)
{
	return 4.0;
}

void hmca_mf_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn);

void hmca_mf_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn);


// Heterogeneous Mean-Field Approximation

#define hmca_hmf_nn_1x1 hmca_mf_nn_1x1

#define hmca_hmf_nn_2x1 hmca_mf_nn_2x1

void hmca_hmf_average (
		const double *y,
		double *average,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights);

void hmca_hmf_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn);

void hmca_hmf_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn);


// Pair Approximation

inline double hmca_pa_nn_1x1 (const int *indices, int n_species_0)
{
	return 3.0;
}

inline double hmca_pa_nn_2x1 (const int *indices, int n_species_0)
{
	return ((indices[0] < n_species_0) == (indices[2] < n_species_0)) ? (2.0) : (4.0);
}

void hmca_pa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn);

void hmca_pa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn);

void hmca_pa_jac_k (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn);


// Half Heterogeneous Pair Approximation

#define hmca_hhpa_nn_1x1 hmca_pa_nn_1x1

#define hmca_hhpa_nn_2x1 hmca_pa_nn_2x1

void hmca_hhpa_average (
		const double *y,
		double *average,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights);

void hmca_hhpa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn);

void hmca_hhpa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn);

void hmca_hhpa_jac_k (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn);


// Machine Learning Moment Closure

void hmca_mlmc_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		mlmc_closure closure, mlmc_closure deriv, void *model);

void hmca_mlmc_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		mlmc_closure closure, mlmc_closure deriv, void *model);

void hmca_mlmc_jac_z (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		mlmc_closure closure, mlmc_closure deriv, void *model);


// Select Pair Approximation

inline double hmca_spa_nn_1x1 (const int *indices, int n_species_0)
{
	return (indices[2] < 0) ? (4.0) : (3.0);
}

inline double hmca_spa_nn_2x1 (const int *indices, int n_species_0)
{
	return (indices[2] < 0) ? (4.0) : (((indices[0] < n_species_0) == (indices[2] < n_species_0)) ? (2.0) : (4.0));
}

void hmca_spa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_pairs, int n_unimol, int n_bimol,
		const int *pairs, const int *reactions, const double *rates, hmca_nn nn);

void hmca_spa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_pairs, int n_unimol, int n_bimol,
		const int *pairs, const int *reactions, const double *rates, hmca_nn nn);

#endif
