#include "hmca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


double closure (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);

double deriv (const int *indices, const double *y, int n_species_0, int n_species_1, void *model);


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
	int n_react = n_unimol+n_bimol;

	double model[] =
	{
		0.590868085, 0, 0.694257733,
		0.25422269, 0.22707202, 0.73582257,
		0.37307309, 0.31246896, 0.70971129,
		0.13338797, 0, 0.90342026,
		-0, 0.18671382, 1.01410331,
		0, 0.04821106, 0.99838797,
		0, 0, 1.00656501,
		0.15427876, 0.15164, 0.87609729,
		0.1186102, 0, 0.91713003,
		-0, 0.02389918, 1.00415731,
		0, 0.04286553, 1.00823407,
		0.499903286, 0, 0.787069522,
		0.10531266, 0.27258614, 0.92357961,
		0.12184353, -0.17975986, 1.07275885,
		0, 0, 0.99488319,
		0, 0, 0.9930671,
		0.68135591, 0.63306404, 0.32875536,
		0.01038036, 0, 0.99843569,
		0.733018271, 0, 0.698173288,
		1.40948279, 1.66883502, -0.52127581,
		0.180848031, 0, 0.896222961,
		0, 0, 0.9940929,
		0.21505666, 0.2086688, 0.78636961,
		0.17155997, 0.18443736, 0.81945774,
		0, 0.01748194, 0.99255996,
		-0, 0.01589982, 1.00384301,
		0, 0.02737791, 1.00381654,
		0.236828019, 0, 0.877571646,
		0.21554713, 0.24413065, 0.74255165,
		0.05503361, -0, 0.95288497,
		0.00703163, 0, 1.00870159,
		0.17912315, 0.21089988, 0.8272617,
		0.366505901, 0, 0.812404352,
		0, 0.01554948, 1.01171213,
		0.16895401, 0.13202032, 0.90025746,
		0, 0.00740506, 1.00347608,
		0.183389056, 0, 0.910350823,
		0.35378287, 0.36487085, 0.68380699,
		0.06620888, 0.11420984, 0.96796872,
		1.85626388, 0, 0.080904137,
		1.11257385, 1.16974333, -0.14325227,
		1.37366021, 0, 0.326232586,
		-0, -0, 1.26725904,
		-0.37556159, -0.18168726, 1.49956227,
		-1.23682736, -0.44928574, 2.28250059,
		-0, 0.22939219, 1.06747606,
		0.12076481, -0, 0.96656824,
		0.05143642, 0.31867052, 0.93841632,
		-0.689032859, -0, 1.55267901,
		-1.11286148, -0.53991111, 2.188268,
		-0, 0.0206679, 1.0855631,
		0.04437183, -0.05297673, 1.0746956,
		0, 0.07675553, 1.00292318,
		1.88096369, 0, -0,
		0.45414086, 0.03363706, 0.62751463,
		0.03504009, -0.08539046, 1.09629108,
		0.58650853, 0.53715067, 0.40503253,
		0, 0, 1.13101261,
		0, -0.06250476, 1.04059936,
		0, 0.33641318, 0.81478217,
		-0.455990453, -0, 1.1806728,
		-0.25204699, -0.3242645, 1.25251916,
		-0, -0, 1.00979961,
		0, 0, 0.9997689,
		0.04057015, 0.05485316, 0.9057319,
		-0, 0.00294996, 0.9916779,
		0, 0, 0.98750815,
		0, 0, 1.00012297,
		0.8830489, 0.79803184, 0.22553724,
		0, 0, 1.0403355,
		-0.11380402, 0, 1.00667254,
		0.000388403033, 0, 1.01033737,
		-0, 0.02386436, 0.9995597,
		0.16505684, 0.13506773, 0.84929317,
		0, 0, 0.99449026,
		1.14748429, 1.14791688, -0.02882161,
		0.25007096, 0.24211745, 0.72798344,
		0.18893313, 0.32322217, 0.7792586,
		0, 0, 1.06499566,
		0, 0, 1.04896272,
		0.55198893, 0.78032318, 0.39538909,
		0, 0, 1.00956466,
		-0.21255444, 0.05985459, 1.0524723,
		-0, -0, 0.88297265,
		-0, -0, 1.00371475,
		0.00919723, -0.09033952, 1.06118126,
		-0.03751221, 0.08545545, 0.92046893,
		0, 0.01825282, 0.97807207,
		0.03219039, -0, 1.01290363,
		-0, 0.04240478, 0.98702282,
		0, 0, 1.01163146,
		-0.01356733, 0.05278499, 0.94420291,
		0, 0.02422268, 1.00144034,
		0, 0.04999917, 1.00263704,
		0, 0.03250425, 1.01448403,
		0, 0, 0.96240552,
		-0, 0.16576962, 0.84809001,
		0.06775046, 0, 0.98363804,
		0.04022324, -0, 0.96015856,
		0, 0, 1.00407474,
		0, 0.01076734, 1.00634207,
		0, 0.07547764, 0.98131921,
		0, 0, 1.0876157,
		0.80263621, 1.19900827, 0.04036355,
		0, 0, 0.84537771,
		1.27367652, 0, 0.38933447,
		0, 0.05663303, 0.92254237,
		1.18188232, 1.08811662, -0.07079617,
		0, 0, 1.14531742,
		0.76341772, 0.954463, 0,
		0.84072711, 0.81326401, 0.30802669,
		0, 0, 1.00062486,
		0.80948528, 0.64197897, 0.34826533,
		-0, 0.08352238, 1.00533686,
		0, 0.31869279, 0.83937558,
		0, 0.01738288, 1.12184869,
		1.66399831, 0, 0.248681264,
		0.98590653, 0.72648186, 0.05454157,
		1.07699766, 1.09797762, -0.36109346,
		1.02386234, 0.49903173, 0.17370906,
		0, 0, 1.24843129,
		0.83671977, 0.85692433, 0.20211261,
		0.43448643, 0.86203777, 0.32764719,
		1.49277609, 0, 0.327362568,
		1.08015374, 1.24281318, -0.28093738,
		2.11972765, 0, -0.215949489,
	};

	double *rates = (double*)malloc(n_react*sizeof(double));
	double *y = (double*)malloc(nn_species*sizeof(double));

	double h = 1e-8;
	double aerr;
	double dfdy[441];
	double dydta[21];
	double dydtb[21];



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
	printf("k =");
	for (i = 0; i < n_react; ++i)
	{
		rates[i] = (double)(rand()%p_kvar+1)/p_kvar;
		printf(" %.1f", rates[i]);
	}
	printf("\n");

	// Random Moments
	printf("y =");
	for (i = 0; i < nn_species; ++i)
	{
		y[i] = (rand()%p_zero > 0) * (double)(rand()%p_yvar+1)/p_yvar;
		printf(" %.2f", y[i]);
	}
	printf("\n");


	// Function and Jacobian
	hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
	hmca_mlmc_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);

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
	printf("J - Jn =\n");
	for (i = 0; i < nn_species; ++i)
	{
		if (y[i] > 0.0)
		{
			y[i] += h;
			hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
			y[i] -= 2.0*h;
			hmca_mlmc_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
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
			hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
			y[i] -= h;
			hmca_mlmc_func(y, dydtb, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
			for(j = 0; j < nn_species; ++j)
			{
				aerr = (dydta[j]-dydtb[j])/h - dfdy[nn_species*j+i];
				printf(" %+.6f", aerr);
			}
		}
		printf("\n");
	}

/*
	// Timing
	t0 = clock();
	for (i = 0; i < 4.3e+4; ++i)
		hmca_mlmc_func(y, dydta, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
	t1 = clock();
	for (i = 0; i < 1.6e+3; ++i)
		hmca_mlmc_jac(y, dfdy, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, closure, deriv, model);
	t2 = clock();
	printf("%f\n", (double)(t1-t0)/CLOCKS_PER_SEC);
	printf("%f\n", (double)(t2-t1)/CLOCKS_PER_SEC);
*/


	free(rates);
	free(y);

	return 0;
}


double closure (const int *indices, const double *y, int n_species_0, int n_species_1, void *model)
{
	double *coef = (double*)model;

	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int ij = hmca_sym_id(indices[0], indices[1], n_species);
	int jk = hmca_sym_id(indices[1], indices[2], n_species);
	int ijk = nn_species*indices[1]+hmca_sym_id(indices[0], indices[2], n_species);

	int i;
	double f0, f1, f2, f3;
	double yi, yj, yk;

	yi = 0.0;
	yj = 0.0;
	yk = 0.0;
	for (i = 0; i < n_species; ++i)
	{
		yi += y[hmca_sym_id(indices[0], i, n_species)];
		yj += y[hmca_sym_id(indices[1], i, n_species)];
		yk += y[hmca_sym_id(indices[2], i, n_species)];
	}

	f0 = yi*yj*yk;
	f1 = y[ij]*yk;
	f2 = yi*y[jk];
	f3 = y[ij]*y[jk]/(yj+(yj == 0.0));

	return (2 + ((indices[0] < n_species_0) != (indices[2] < n_species_0)) * 2) * (
			(1.0-coef[3*ijk+0]-coef[3*ijk+1]-coef[3*ijk+2]) * f0
			+ coef[3*ijk+0] * f1
			+ coef[3*ijk+1] * f2
			+ coef[3*ijk+2] * f3
			);
}


double deriv (const int *indices, const double *y, int n_species_0, int n_species_1, void *model)
{
	double *coef = (double*)model;

	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int ij = hmca_sym_id(indices[0], indices[1], n_species);
	int jk = hmca_sym_id(indices[1], indices[2], n_species);
	int ijk = nn_species*indices[1]+hmca_sym_id(indices[0], indices[2], n_species);
	int ab = hmca_sym_id(indices[3], indices[4], n_species);

	int i;
	double df0, df1, df2, df3;
	double yi, yj, yk;

	yi = 0.0;
	yj = 0.0;
	yk = 0.0;
	for (i = 0; i < n_species; ++i)
	{
		yi += y[hmca_sym_id(indices[0], i, n_species)];
		yj += y[hmca_sym_id(indices[1], i, n_species)];
		yk += y[hmca_sym_id(indices[2], i, n_species)];
	}

	df0 = (indices[0] == indices[3] || indices[0] == indices[4]) * yj*yk
		+ (indices[1] == indices[3] || indices[1] == indices[4]) * yi*yk
		+ (indices[2] == indices[3] || indices[2] == indices[4]) * yi*yj;
	df1 = (ab == ij) * yk + (indices[2] == indices[3] || indices[2] == indices[4]) * y[ij];
	df2 = (ab == jk) * yi + (indices[0] == indices[3] || indices[0] == indices[4]) * y[jk];
	df3 = (ab == ij) * (y[jk]+(yj == 0.0 && ab == jk))/(yj+(yj == 0.0))
		+ (ab == jk) * (y[ij]+(yj == 0.0 && ab == ij))/(yj+(yj == 0.0))
		- (indices[1] == indices[3] || indices[1] == indices[4]) * (y[ij]*y[jk]+(yj == 0.0 && ab == ij && ab == jk))/(yj*yj+(yj == 0.0));

	return (2 + ((indices[0] < n_species_0) != (indices[2] < n_species_0)) * 2) * (
			(1.0-coef[3*ijk+0]-coef[3*ijk+1]-coef[3*ijk+2]) * df0
			+ coef[3*ijk+0] * df1
			+ coef[3*ijk+1] * df2
			+ coef[3*ijk+2] * df3
			);
}
