#include "hmca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


typedef struct
{
	int idz;
	double dz;
}
data_Jz;


double closure (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);

double deriv (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);

double closure_dz (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);


int main (int argc, char* argv[])
{
	clock_t t0, t1, t2;
	int i, j, k, l;

	// Reaction Network
	enum species {Vc, Cc, Oc, Vb, Cb, Ob};
	int n_species_0 = 3;
	int n_species_1 = 3;
	int n_unimol = 4;
	int n_bimol = 18;
	int reactions[] =
	{
		Vc,Cc, Vb,Cb, Cc,Vc, Cb,Vb,
		Vc,Vc,Oc,Oc, Vb,Vb,Ob,Ob, Vb,Vc,Ob,Oc, Oc,Oc,Vc,Vc, Ob,Ob,Vb,Vb, Ob,Oc,Vb,Vc,
		Cc,Vc,Vc,Cc, Cb,Vb,Vb,Cb, Cc,Vb,Vc,Cb, Cb,Vc,Vb,Cc,
		Oc,Vc,Vc,Oc, Ob,Vb,Vb,Ob, Oc,Vb,Vc,Ob, Ob,Vc,Vb,Oc,
		Cc,Oc,Vc,Vc, Cb,Ob,Vb,Vb, Cc,Ob,Vc,Vb, Cb,Oc,Vb,Vc,
	};
/*	double vals[] =
	{
		7.2e+8, 7.2e+8, 9.2e+6, 2.8e+4,
		9.7e+7, 9.7e+7, 9.7e+7, 2.8e+1, 4.1e-21, 3.4e-10,
		6.6e-2, 1.1e+8, 1.5e+2, 0.5,
		0.5, 1.6e+7, 4.9e+4, 6.0e-7,
		1.7e+5, 1.6, 5.2e+2, 1.2e+6,
	};
	double concs[] =
	{
		0.5, 0.0, 0.0, 0.5, 0.0, 0.0,
	};	*/

	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int nnn_species = n_species*nn_species;
	int n_react = n_unimol+n_bimol;

	double *rates = (double*)malloc(n_react*sizeof(double));
	double *y = (double*)malloc(nn_species*sizeof(double));

	double h = 1e-8;
	double aerr;
	double *dfdy = (double*)malloc(nn_species*nn_species*sizeof(double));
	double *dydta = (double*)malloc(nn_species*sizeof(double));
	double *dydtb = (double*)malloc(nn_species*sizeof(double));
	double *dfdz = (double*)malloc(nn_species*nnn_species*sizeof(double));


	// Random Number Generator
	unsigned long seed;
	int p_kvar = 10;
	int p_zero = 3;
	int p_yvar = 100;
	if (argc > 1)
		sscanf(argv[1], "%ld", &seed);
	else
		seed = (unsigned long)time(NULL);
	srand(seed);
	printf("seed = %ld\n", seed);

	// Random Rate Constants
	for (i = 0; i < n_react; ++i)
		rates[i] = (double)(rand()%p_kvar+1)/p_kvar;

	printf("k =");
	for (i = 0; i < n_react; ++i)
		printf(" %.1f", rates[i]);
	printf("\n");

	// Random Moments
	for (i = 0; i < nn_species; ++i)
		y[i] = (rand()%p_zero > 0) * (double)(rand()%p_yvar+1)/p_yvar;

	printf("y =");
	for (i = 0; i < nn_species; ++i)
		printf(" %.2f", y[i]);
	printf("\n");


	// Function and Jacobian
	hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
	hmca_mlmc_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);

	printf("f =");
	for (i = 0; i < nn_species; ++i)
		printf(" %+.6f", dydta[i]);
	printf("\n");

	printf("J =\n");
	for (i = 0; i < nn_species; ++i)
	{
		for(j = 0; j < nn_species; ++j)
			printf(" %+.6f", dfdy[nn_species*i+j]);
		printf("\n");
	}

	// Compare to Numerical Jacobian
	printf("J_n - J =\n");
	for (i = 0; i < nn_species; ++i)
	{
		if (y[i] > 0.0)
		{
			y[i] += h;
			hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
			y[i] -= 2.0*h;
			hmca_mlmc_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
			y[i] += h;
			for(j = 0; j < nn_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdy[nn_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		else
		{
			y[i] += h;
			hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
			y[i] -= h;
			hmca_mlmc_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
			for(j = 0; j < nn_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/h - dfdy[nn_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		printf("\n");
	}

	// Jacobian w.r.t. Closure
	hmca_mlmc_dfdz(y, dfdz, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);

	printf("dfdz_n - dfdz =\n");
	for (i = 0; i < nnn_species; ++i)
	{
		data_Jz forward = {i, +h};
		hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure_dz, deriv, &forward);
		data_Jz backward = {i, -h};
		hmca_mlmc_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure_dz, deriv, &backward);
		for(j = 0; j < nn_species; ++j)
		{
			aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdz[nnn_species*j+i];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}

/*	// Timing
	t0 = clock();
	for (i = 0; i < 9.0e+4; ++i)
		hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
	t1 = clock();
	for (i = 0; i < 3.2e+3; ++i)
		hmca_mlmc_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mlmc_nn_2x1, closure, deriv, NULL);
	t2 = clock();
	printf("%f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	printf("%f\n", (double)(t2-t1)/CLOCKS_PER_SEC);	*/


	free(rates);
	free(y);

	free(dfdy);
	free(dydta);
	free(dydtb);
	free(dfdz);

	return 0;
}


double closure (const int *indices, const double *y, int n_species_0, int n_species_1, void *model)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int ij = hmca_sym_id(indices[0], indices[1], n_species);
	int jk = hmca_sym_id(indices[1], indices[2], n_species);
	int ijk = nn_species*indices[1]+hmca_sym_id(indices[0], indices[2], n_species);

	int a;
	double yj;

	yj = 0.0;
	for (a = 0; a < n_species; ++a)
		yj += y[hmca_sym_id(indices[1], a, n_species)];

	return y[ij]*y[jk]/(yj+(yj == 0.0));
}


double deriv (const int *indices, const double *y, int n_species_0, int n_species_1, void *model)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int ij = hmca_sym_id(indices[0], indices[1], n_species);
	int jk = hmca_sym_id(indices[1], indices[2], n_species);
	int ijk = nn_species*indices[1]+hmca_sym_id(indices[0], indices[2], n_species);
	int lm = hmca_sym_id(indices[3], indices[4], n_species);

	int a;
	double yj;

	yj = 0.0;
	for (a = 0; a < n_species; ++a)
		yj += y[hmca_sym_id(indices[1], a, n_species)];

	return (lm == ij) * (y[jk]+(yj == 0.0 && lm == jk))/(yj+(yj == 0.0))
		+ (lm == jk) * (y[ij]+(yj == 0.0 && lm == ij))/(yj+(yj == 0.0))
		- (indices[1] == indices[3] || indices[1] == indices[4]) * (y[ij]*y[jk]+(yj == 0.0 && lm == ij && lm == jk))/(yj*yj+(yj == 0.0));
}


double closure_dz (const int *indices, const double *y, int n_species_0, int n_species_1, void *model)
{
	data_Jz *p = (data_Jz*)model;

	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int ij = hmca_sym_id(indices[0], indices[1], n_species);
	int jk = hmca_sym_id(indices[1], indices[2], n_species);
	int ijk = nn_species*indices[1]+hmca_sym_id(indices[0], indices[2], n_species);

	int a;
	double yj;

	yj = 0.0;
	for (a = 0; a < n_species; ++a)
		yj += y[hmca_sym_id(indices[1], a, n_species)];

	return y[ij]*y[jk]/(yj+(yj == 0.0)) + (ijk == p->idz) * p->dz;
}
