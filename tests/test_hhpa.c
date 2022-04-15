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
	int mesh = 3;
	double bound = 1.0;

	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *logk0 = (double*)malloc(n_react*sizeof(double));
	double *dlogk = (double*)malloc(n_react*sizeof(double));
	double *rates = (double*)malloc(mesh*n_react*sizeof(double));
	double *weights = (double*)malloc(mesh*sizeof(double));
	double *y = (double*)malloc(mesh*n_species*n_species*sizeof(double));

	double h = 1e-8;
	double aerr, rerr;
	double *dfdy = (double*)malloc(mesh*n_species*n_species*mesh*n_species*n_species*sizeof(double));
	double *dydta = (double*)malloc(mesh*n_species*n_species*sizeof(double));
	double *dydtb = (double*)malloc(mesh*n_species*n_species*sizeof(double));
	double *dfdk = (double*)malloc(mesh*n_species*n_species*mesh*n_react*sizeof(double));
	double *dkdlogk0 = (double*)malloc(mesh*n_react*sizeof(double));
	double *dkddlogk = (double*)malloc(mesh*n_react*sizeof(double));
	double *ka = (double*)malloc(mesh*n_react*sizeof(double));
	double *kb = (double*)malloc(mesh*n_react*sizeof(double));
	double *wt = (double*)malloc(mesh*sizeof(double));
	double *dfdlogk0 = (double*)malloc(mesh*n_species*n_species*n_react*sizeof(double));
	double *dfddlogk = (double*)malloc(mesh*n_species*n_species*n_react*sizeof(double));


	// Random Number Generator
	unsigned long seed;
	int p_kvar = 10;
	int p_zero = 5;
	int p_yvar = 100;
	if (argc > 1)
		sscanf(argv[1], "%ld", &seed);
	else
		seed = (unsigned long)time(NULL);
	srand(seed);
	printf("seed = %ld\n", seed);

	// Random Rate Constants
	for (i = 0; i < mesh; ++i) for (j = 0; j < n_react; ++j)
		rates[n_react*i+j] = (double)(rand()%p_kvar+1)/p_kvar;

	printf("k =\n");
	for (i = 0; i < mesh; ++i)
	{
		for (j = 0; j < n_react; ++j)
			printf(" %+.2f", rates[n_react*i+j]);
		printf("\n");
	}

	// Random Weights
	for (i = 0; i < mesh; ++i)
		weights[i] = (double)(rand()%p_kvar+1)/p_kvar;

	printf("w =\n");
	for (i = 0; i < mesh; ++i)
		printf(" %+.2f", weights[i]);
	printf("\n");

/*	// Random Rate Constants (Known Distribution)
	for (i = 0; i < n_react; ++i)
	{
		logk0[i] = log((double)(rand()%p_kvar+1)/p_kvar);
		dlogk[i] = rand()%p_zero;
	}
	hmca_lognorm_set(logk0, dlogk, rates, weights, n_unimol, n_bimol, mesh, bound);

	printf("k =\n");
	for (i = 0; i < mesh; ++i)
	{
		for (j = 0; j < n_react; ++j)
			printf(" %+.2f", rates[n_react*i+j]);
		printf("\n");
	}

	printf("w =\n");
	for (i = 0; i < mesh; ++i)
		printf(" %+.2f", weights[i]);
	printf("\n");	*/

	// Random Moments
	for (i = 0; i < n_species; ++i) for (j = i; j < n_species; ++j)
		y[n_species*i+j] = (rand()%p_zero > 0) * (double)(rand()%p_yvar+1)/p_yvar;
	for (i = 0; i < n_species; ++i) for (j = 0; j < i; ++j)
		y[n_species*i+j] = (y[n_species*j+i] > 0.0) * (double)(rand()%p_yvar+1)/p_yvar;
	for (i = 1; i < mesh; ++i) for (j = 0; j < n_species*n_species; ++j)
		y[n_species*n_species*i+j] = (y[j] > 0.0) * (double)(rand()%p_yvar+1)/p_yvar;

	printf("y =\n");
	for (i = 0; i < mesh; ++i)
	{
		for (j = 0; j < n_species*n_species; ++j)
			printf(" %+.2f", y[n_species*n_species*i+j]);
		printf("\n");
	}

	// Function and Jacobian
	hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
	hmca_hhpa_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);

	printf("f =\n");
	for (i = 0; i < mesh; ++i)
	{
		for (j = 0; j < n_species*n_species; ++j)
			printf(" %+.6f", dydta[n_species*n_species*i+j]);
		printf("\n");
	}

	printf("J =\n");
	for (i = 0; i < mesh*n_species*n_species; ++i)
	{
		for(j = 0; j < mesh*n_species*n_species; ++j)
			printf(" %+.6f", dfdy[mesh*n_species*n_species*i+j]);
		printf("\n");
	}

/*	// Compare to Numerical Jacobian
	printf("J_n - J =\n");
	for (i = 0; i < mesh*n_species*n_species; ++i)
	{

		if (y[i] > 0.0)
		{
			y[i] += h;
			hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
			y[i] -= 2.0*h;
			hmca_hhpa_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
			y[i] += h;
			for(j = 0; j < mesh*n_species*n_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdy[mesh*n_species*n_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		else
		{
			y[i] += h;
			hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
			y[i] -= h;
			hmca_hhpa_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
			for(j = 0; j < mesh*n_species*n_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/h - dfdy[mesh*n_species*n_species*j+i];
				printf(" %+.6f", aerr);
			}

		}
		printf("\n");
	}	*/

/*	// Derivative w.r.t. Rate Constants
	hmca_hhpa_dfdk(y, dfdk, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);

	printf("dfdk_n - dfdk =\n");
	for (i = 0; i < mesh*n_react; ++i)
	{
		rates[i] += h;
		hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
		rates[i] -= 2.0*h;
		hmca_hhpa_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_hhpa_nn_2x1);
		rates[i] += h;
		for(j = 0; j < mesh*n_species*n_species; ++j)
		{
			aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdk[mesh*n_react*j+i];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}

	hmca_lognorm_dkdx(logk0, dlogk, dkdlogk0, dkddlogk, n_unimol, n_bimol, mesh, bound);

	printf("dkdlogk0 =\n");
	for (i = 0; i < mesh*n_react; ++i)
		printf(" %+.6f", dkdlogk0[i]);
	printf("\n");

	printf("dkddlogk =\n");
	for (i = 0; i < mesh*n_react; ++i)
		printf(" %+.6f", dkddlogk[i]);
	printf("\n");

	printf("dkdlogk0_n - dkdlogk0 =\n");
	for (i = 0; i < n_react; ++i)
	{
		logk0[i] += h;
		hmca_lognorm_set(logk0, dlogk, ka, wt, n_unimol, n_bimol, mesh, bound);
		logk0[i] -= 2.0*h;
		hmca_lognorm_set(logk0, dlogk, kb, wt, n_unimol, n_bimol, mesh, bound);
		logk0[i] += h;
		for(j = 0; j < mesh*n_react; ++j)
		{
			aerr = (ka[j]-kb[j])/(2.0*h) - (j%n_react == i) * dkdlogk0[j];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}
	
	printf("dkddlogk_n - dkddlogk =\n");
	for (i = 0; i < n_react; ++i)
	{
		dlogk[i] += h;
		hmca_lognorm_set(logk0, dlogk, ka, wt, n_unimol, n_bimol, mesh, bound);
		dlogk[i] -= 2.0*h;
		hmca_lognorm_set(logk0, dlogk, kb, wt, n_unimol, n_bimol, mesh, bound);
		dlogk[i] += h;
		for(j = 0; j < mesh*n_react; ++j)
		{
			aerr = (ka[j]-kb[j])/(2.0*h) - (j%n_react == i) * dkddlogk[j];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}

	memset(dfdlogk0, 0, mesh*n_species*n_species*n_react*sizeof(double));
	for (i = 0; i < mesh*n_species*n_species; ++i) for (j = 0; j < mesh*n_react; ++j) for (k = 0; k < n_react; ++k)
		dfdlogk0[n_react*i+k] += (j%n_react == k%n_react) * dfdk[mesh*n_react*i+j] * dkdlogk0[j];

	memset(dfddlogk, 0, mesh*n_species*n_species*n_react*sizeof(double));
	for (i = 0; i < mesh*n_species*n_species; ++i) for (j = 0; j < mesh*n_react; ++j) for (k = 0; k < n_react; ++k)
		dfddlogk[n_react*i+k] += (j%n_react == k%n_react) * dfdk[mesh*n_react*i+j] * dkddlogk[j];

	printf("dfdlogk0_n - dfdlogk0 =\n");
	for (i = 0; i < n_react; ++i)
	{
		logk0[i] += h;
		hmca_lognorm_set(logk0, dlogk, ka, wt, n_unimol, n_bimol, mesh, bound);
		hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, ka, weights, hmca_hhpa_nn_2x1);
		logk0[i] -= 2.0*h;
		hmca_lognorm_set(logk0, dlogk, kb, wt, n_unimol, n_bimol, mesh, bound);
		hmca_hhpa_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, kb, weights, hmca_hhpa_nn_2x1);
		logk0[i] += h;
		for(j = 0; j < mesh*n_species*n_species; ++j)
		{
			aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfdlogk0[n_react*j+i];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}

	printf("dfddlogk_n - dfddlogk =\n");
	for (i = 0; i < n_react; ++i)
	{
		dlogk[i] += h;
		hmca_lognorm_set(logk0, dlogk, ka, wt, n_unimol, n_bimol, mesh, bound);
		hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, ka, weights, hmca_hhpa_nn_2x1);
		dlogk[i] -= 2.0*h;
		hmca_lognorm_set(logk0, dlogk, kb, wt, n_unimol, n_bimol, mesh, bound);
		hmca_hhpa_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, kb, weights, hmca_hhpa_nn_2x1);
		dlogk[i] += h;
		for(j = 0; j < mesh*n_species*n_species; ++j)
		{
			aerr = (dydta[j]-dydtb[j])/(2.0*h) - dfddlogk[n_react*j+i];
			printf(" %+.6f", aerr);
		}
		printf("\n");
	}	*/

/*	// Timing
	t0 = clock();
	for (i = 0; i < 5.0e+4; ++i)
		hmca_hhpa_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_mf_nn_2x1);
	t1 = clock();
	for (i = 0; i < 2.5e+3; ++i)
		hmca_hhpa_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, hmca_mf_nn_2x1);
	t2 = clock();
	printf("%f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	printf("%f\n", (double)(t2-t1)/CLOCKS_PER_SEC);	*/


	free(logk0);
	free(dlogk);
	free(rates);
	free(weights);
	free(y);

	free(dfdy);
	free(dydta);
	free(dydtb);
	free(dfdk);
	free(dkdlogk0);
	free(dkddlogk);
	free(ka);
	free(kb);
	free(wt);
	free(dfdlogk0);
	free(dfddlogk);

	return 0;
}
