#include "hmca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


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
	int n_react = n_unimol+n_bimol;

	double *rates = (double*)malloc(n_react*sizeof(double));
	double *y = (double*)malloc(n_species*sizeof(double));

	double h = 1e-8;
	double aerr;
	double *dfdy = (double*)malloc(n_species*n_species*sizeof(double));
	double *dydta = (double*)malloc(n_species*sizeof(double));
	double *dydtb = (double*)malloc(n_species*sizeof(double));
	double *dfdk = (double*)malloc(n_species*n_react*sizeof(double));


	// Randum Number Generator
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
	for (i = 0; i < n_species; ++i)
		y[i] = (rand()%p_zero > 0) * (double)(rand()%p_yvar+1)/p_yvar;

	printf("y =");
	for (i = 0; i < n_species; ++i)
		printf(" %.2f", y[i]);
	printf("\n");


	// Function and Jacobian
	hmca_mf_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
	hmca_mf_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);

	printf("f =");
	for (i = 0; i < n_species; ++i)
		printf(" %+.6f", dydta[i]);
	printf("\n");

	printf("J =\n");
	for (i = 0; i < n_species; ++i)
	{
		for(j = 0; j < n_species; ++j)
			printf(" %+.6f", dfdy[n_species*i+j]);
		printf("\n");
	}

	// Compare to Numerical Jacobian
	printf("J_n - J =\n");
	for (i = 0; i < n_species; ++i)
	{
		if (y[i] > 0.0)
		{
			y[i] += h;
			hmca_mf_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
			y[i] -= 2.0*h;
			hmca_mf_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
			y[i] += h;
			for(j = 0; j < n_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdy[n_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		else
		{
			y[i] += h;
			hmca_mf_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
			y[i] -= h;
			hmca_mf_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
			for(j = 0; j < n_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/h - dfdy[n_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		printf("\n");
	}

	// Derivative w.r.t. Rate Constants
	hmca_mf_dfdk(y, dfdk, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);

	printf("dfdk_n - dfdk =\n");
	for (i = 0; i < n_react; ++i)
	{
		rates[i] += h;
		hmca_mf_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
		rates[i] -= 2.0*h;
		hmca_mf_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
		rates[i] += h;
		for(j = 0; j < n_species; ++j)
		{
			aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdk[n_react*j+i];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}


/*	// Timing
	t0 = clock();
	for (i = 0; i < 4.2e+6; ++i)
		hmca_mf_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
	t1 = clock();
	for (i = 0; i < 2.2e+6; ++i)
		hmca_mf_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca_mf_nn_2x1);
	t2 = clock();
	printf("%f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	printf("%f\n", (double)(t2-t1)/CLOCKS_PER_SEC);	*/


	free(rates);
	free(y);

	free(dfdy);
	free(dydta);
	free(dydtb);
	free(dfdk);

	return 0;
}
