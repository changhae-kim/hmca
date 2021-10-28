#include "hmca.h"


// Shared Functions

extern inline int hmca_sym_id (int i, int j, int n);

void hmca_lognorm_set (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;
	double norm;

	/*
	// midpoint rule
	dx = 2.0 * bound / (mesh + (mesh == 0));
	for (i = 0; i < mesh; ++i)
	{
		x = (i+0.5)*dx-bound;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]+x*dlogk[j]);
		weights[i] = exp(-0.5*x*x)*dx;
	}
	*/

	// trapezoid rule
	dx = 2.0 * bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx-bound;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]+x*dlogk[j]);
		weights[i] = (1.0 + (i > 0 && i+1 < mesh) * 1.0) * exp(-0.5*x*x) * dx;
	}

	norm = 0.0;
	for (i = 0; i < mesh; ++i)
		norm += weights[i];
	for (i = 0; i < mesh; ++i)
		weights[i] /= norm;

	return;
}

void hmca_lognorm_dkdx (
		const double *logk0, const double *dlogk,
		double *dkdlogk0, double *dkddlogk,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;

	/*
	// midpoint rule
	dx = 2.0 * bound / (mesh + (mesh == 0));
	for (i = 0; i < mesh; ++i)
	{
		x = (i+0.5)*dx-bound;
		for (j = 0; j < n_react; ++j)
		{
			dkdlogk0[n_react*i+j] = exp(logk0[j]+x*dlogk[j]);
			dkddlogk[n_react*i+j] = x*exp(logk0[j]+x*dlogk[j]);
		}
	}
	*/

	//  trapezoid rule
	dx = 2.0 * bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx-bound;
		for (j = 0; j < n_react; ++j)
		{
			dkdlogk0[n_react*i+j] = exp(logk0[j]+x*dlogk[j]);
			dkddlogk[n_react*i+j] = x*exp(logk0[j]+x*dlogk[j]);
		}
	}

	return;
}

void hmca_logexp_set (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;
	double norm;

	/*
	// left rule
	dx = bound / (mesh + (mesh == 0));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]-x*dlogk[j]);
		weights[i] = exp(-x)*dx;
	}
	*/

	// Simpson's rule
	dx = bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]-x*dlogk[j]);
		weights[i] = (1.0 + (i > 0 && i+1 < mesh) * (1.0 + (i%2 == 1) * 2.0)) * exp(-x) * dx;
	}

	norm = 0.0;
	for (i = 0; i < mesh; ++i)
		norm += weights[i];
	for (i = 0; i < mesh; ++i)
		weights[i] /= norm;

	return;
}

void hmca_logexp_dkdx (
		const double *logk0, const double *dlogk,
		double *dkdlogk0, double *dkddlogk,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;
	double norm;

	/*
	// left rule
	dx = bound / (mesh + (mesh == 0));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx;
		for (j = 0; j < n_react; ++j)
		{
			dkdlogk0[n_react*i+j] = exp(logk0[j]-x*dlogk[j]);
			dkddlogk[n_react*i+j] = -x*exp(logk0[j]-x*dlogk[j]);
		}
	}
	*/

	// Simpson's rule
	dx = bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx;
		for (j = 0; j < n_react; ++j)
		{
			dkdlogk0[n_react*i+j] = exp(logk0[j]-x*dlogk[j]);
			dkddlogk[n_react*i+j] = -x*exp(logk0[j]-x*dlogk[j]);
		}
	}

	return;
}

void hmca_logsech_set (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;

	// midpoint rule
//	dx = 2.0 * bound / (mesh + (mesh == 0));
//	for (i = 0; i < mesh; ++i)
//	{
//		x = (i+0.5)*dx-bound;
//		for (j = 0; j < n_react; ++j)
//			rates[n_react*i+j] = exp(logk[j]+x*dlogk[j]);
//		weights[i] = (2.0/M_PI) / (exp(+x)+exp(-x)) * dx;
//	}

	// trapezoid rule
	dx = 2.0 * bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx-bound;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]+x*dlogk[j]);
		weights[i] = (1 + (i > 0 && i+1 < mesh)) / 2.0 * (2.0/M_PI) / (exp(+x)+exp(-x)) * dx;
	}

	dx = 0.0;
	for (i = 0; i < mesh; ++i)
		dx += weights[i];
	for (i = 0; i < mesh; ++i)
		weights[i] /= dx;

	return;
}

void hmca_logpoisson2_set (
		const double *logk0, const double *dlogk,
		double *rates, double *weights,
		int n_unimol, int n_bimol, int mesh, double bound
		)
{
	int n_react = n_unimol+n_bimol;

	int i, j;
	double x, dx;

	// midpoint rule
//	dx = bound / (mesh + (mesh == 0));
//	for (i = 0; i < mesh; ++i)
//	{
//		x = (i+0.5)*dx;
//		for (j = 0; j < n_react; ++j)
//			rates[n_react*i+j] = exp(logk[j]-x*dlogk[j]);
//		weights[i] = x * x * exp(-x) * dx;
//	}

	// Simpson's rule
	dx = bound / (mesh-1 + (mesh == 1));
	for (i = 0; i < mesh; ++i)
	{
		x = i*dx;
		for (j = 0; j < n_react; ++j)
			rates[n_react*i+j] = exp(logk0[j]-x*dlogk[j]);
		weights[i] = (1 + (i > 0 && i+1 < mesh) * (1 + (i%2 == 1) * 2)) / 3.0 * x * x * exp(-x) * dx;
	}

	dx = 0.0;
	for (i = 0; i < mesh; ++i)
		dx += weights[i];
	for (i = 0; i < mesh; ++i)
		weights[i] /= dx;

	return;
}


// Mean-Field Approximation

extern inline double hmca_mf_nn_1x1 (const int *indices, int n_species_0);

extern inline double hmca_mf_nn_2x1 (const int *indices, int n_species_0);

void hmca_mf_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int i;
	int r0, r1, p0, p1;
	int indices[2];
	double ky;

	memset(dydt, 0, (n_species_0+n_species_1)*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		ky = rates[i] * y[r0];
		dydt[r0] -= ky;
		dydt[p0] += ky;
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		indices[0] = r0;
		indices[1] = r1;
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		dydt[r0] -= ky;
		dydt[p0] += ky;

		indices[0] = r1;
		indices[1] = r0;
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		dydt[r1] -= ky;
		dydt[p1] += ky;
	}

	return;
}

void hmca_mf_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;

	int i;
	int r0, r1, p0, p1;
	int indices[2];
	double ky;

	memset(dfdy, 0, n_species*n_species*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		// ky = rates[i] * y[r0];
		dfdy[n_species*r0+r0] -= rates[i];
		dfdy[n_species*p0+r0] += rates[i];
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		indices[0] = r0;
		indices[1] = r1;
		// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r1];
		dfdy[n_species*r0+r0] -= ky;
		dfdy[n_species*p0+r0] += ky;
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0];
		dfdy[n_species*r0+r1] -= ky;
		dfdy[n_species*p0+r1] += ky;

		indices[0] = r1;
		indices[1] = r0;
		// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r1];
		dfdy[n_species*r1+r0] -= ky;
		dfdy[n_species*p1+r0] += ky;
		ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0];
		dfdy[n_species*r1+r1] -= ky;
		dfdy[n_species*p1+r1] += ky;
	}

	return;
}

void hmca_mf_dfdk (
		const double *y,
		double *dfdk,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_react = n_unimol+n_bimol;

	int i;
	int r0, r1, p0, p1;
	int indices[2];
	double ky;

	memset(dfdk, 0, (n_species_0+n_species_1)*n_react*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		// ky = rates[i] * y[r0];
		dfdk[n_react*r0+i] -= y[r0];
		dfdk[n_react*p0+i] += y[r0];
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		indices[0] = r0;
		indices[1] = r1;
		// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		ky = (*nn) (indices, n_species_0) * y[r0] * y[r1];
		dfdk[n_react*r0+n_unimol+i] -= ky;
		dfdk[n_react*p0+n_unimol+i] += ky;

		indices[0] = r1;
		indices[1] = r0;
		// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
		ky = (*nn) (indices, n_species_0) * y[r0] * y[r1];
		dfdk[n_react*r1+n_unimol+i] -= ky;
		dfdk[n_react*p1+n_unimol+i] += ky;
	}

	return;
}


// Heterogeneous Mean-Field Approximation

void hmca_hmf_average (
		const double *y,
		double *_y, double *_ky,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	memset(_y, 0, n_species*sizeof(double));
	memset(_ky, 0, n_bimol*sizeof(double));

	for (int j = 0; j < mesh; ++j)
	{
		const double *kj = rates+n_react*j;
		const double *yj = y+n_species*j;

		int k;
		int r0;

		for (k = 0; k < n_species; ++k)
		{
			_y[k] += yj[k] * weights[j];
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			_ky[k] += kj[n_unimol+k] * yj[r0] * weights[j];
		}
	}

	return;
}

void hmca_hmf_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*sizeof(double));
	double *_ky = (double*)malloc(n_bimol*sizeof(double));
	hmca_hmf_average(
			y,
			_y, _ky,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights);

	for (int i = 0; i < mesh; ++i)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*i;
		double *fi = dydt+n_species*i;

		int k;
		int r0, r1, p0, p1;
		int indices[2];
		double ky;

		memset(fi, 0, n_species*sizeof(double));

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			p0 = reactions[2*k+1];
			ky = ki[k] * yi[r0];
			fi[r0] -= ky;
			fi[p0] += ky;
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			indices[0] = r0;
			indices[1] = r1;
			ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * _y[r1];
			fi[r0] -= ky;
			fi[p0] += ky;

			indices[0] = r1;
			indices[1] = r0;
			ky = (*nn) (indices, n_species_0) * _ky[k] * yi[r1];
			fi[r1] -= ky;
			fi[p1] += ky;
		}
	}

	free(_y);
	free(_ky);

	return;
}

void hmca_hmf_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*sizeof(double));
	double *_ky = (double*)malloc(n_bimol*sizeof(double));
	hmca_hmf_average(
			y,
			_y, _ky,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights);

	for (int i = 0; i < mesh; ++i) for (int j = 0; j < mesh; ++j)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*i;
		const double *kj = rates+n_react*j;
		double **dfij = (double**)malloc(n_species*sizeof(double*));

		int k;
		int r0, r1, p0, p1;
		int indices[2];
		double ky;

		for (k = 0; k < n_species; ++k)
		{
			dfij[k] = dfdy + n_species*mesh*(n_species*i+k) + n_species*j;
			memset(dfij[k], 0, n_species*sizeof(double));
		}

		if (i == j)
		{
			for (k = 0; k < n_unimol; ++k)
			{
				r0 = reactions[2*k+0];
				p0 = reactions[2*k+1];
				// ky = ki[k] * yi[r0];
				dfij[r0][r0] -= ki[k];
				dfij[p0][r0] += ki[k];
			}

			for (k = 0; k < n_bimol; ++k)
			{
				r0 = reactions[2*n_unimol+4*k+0];
				r1 = reactions[2*n_unimol+4*k+1];
				p0 = reactions[2*n_unimol+4*k+2];
				p1 = reactions[2*n_unimol+4*k+3];

				indices[0] = r0;
				indices[1] = r1;
				// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * _y[r1];
				ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * _y[r1];
				dfij[r0][r0] -= ky;
				dfij[p0][r0] += ky;

				indices[0] = r1;
				indices[1] = r0;
				// ky = (*nn) (indices, n_species_0) * _ky[k] * yi[r1];
				ky = (*nn) (indices, n_species_0) * _ky[k];
				dfij[r1][r1] -= ky;
				dfij[p1][r1] += ky;
			}
		} // end if (i == j)

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			indices[0] = r0;
			indices[1] = r1;
			// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * _y[r1];
			ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * weights[j];
			dfij[r0][r1] -= ky;
			dfij[p0][r1] += ky;

			indices[0] = r1;
			indices[1] = r0;
			// ky = (*nn) (indices, n_species_0) * _ky[k] * yi[r1];
			ky = (*nn) (indices, n_species_0) * kj[n_unimol+k] * weights[j] * yi[r1];
			dfij[r1][r0] -= ky;
			dfij[p1][r0] += ky;
		}

		free(dfij);
	} // end for (i) for (j)

	free(_y);
	free(_ky);

	return;
}

void hmca_hmf_dfdk (
		const double *y,
		double *dfdk,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*sizeof(double));
	double *_ky = (double*)malloc(n_bimol*sizeof(double));
	hmca_hmf_average(
			y,
			_y, _ky,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights
			);

	for (int i = 0; i < mesh; ++i) for (int j = 0; j < mesh; ++j)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*i;
		const double *kj = rates+n_react*j;
		const double *yj = y+n_species*j;
		double **dfij = (double**)malloc(n_species*sizeof(double*));

		int k, l, m;
		int r0, r1, p0, p1;
		int indices[3];
		double ky;

		for (k = 0; k < n_species; ++k)
		{
			dfij[k] = dfdk + mesh*n_react*(n_species*i+k) + n_react*j;
			memset(dfij[k], 0, n_react*sizeof(double));
		}

		if (i == j)
		{
			for (k = 0; k < n_unimol; ++k)
			{
				r0 = reactions[2*k+0];
				p0 = reactions[2*k+1];
				// ky = ki[k] * yi[r0];
				dfij[r0][k] -= yi[r0];
				dfij[p0][k] += yi[r0];
			}

			for (k = 0; k < n_bimol; ++k)
			{
				r0 = reactions[2*n_unimol+4*k+0];
				r1 = reactions[2*n_unimol+4*k+1];
				p0 = reactions[2*n_unimol+4*k+2];
				p1 = reactions[2*n_unimol+4*k+3];

				indices[0] = r0;
				indices[1] = r1;
				// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * _y[r1];
				ky = (*nn) (indices, n_species_0) * yi[r0] * _y[r1];
				dfij[r0][n_unimol+k] -= ky;
				dfij[p0][n_unimol+k] += ky;

				// indices[0] = r1;
				// indices[1] = r0;
				// ky = (*nn) (indices, n_species_0) * _ky[k] * yi[r1];
			}
		} // end if (i == j)

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			// indices[0] = r0;
			// indices[1] = r1;
			// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[r0] * _y[r1];

			indices[0] = r1;
			indices[1] = r0;
			// ky = (*nn) (indices, n_species_0) * _ky[k] * yi[r1];
			ky = (*nn) (indices, n_species_0) * yj[r0] * weights[j] * yi[r1];
			dfij[r1][n_unimol+k] -= ky;
			dfij[p1][n_unimol+k] += ky;
		}
	
		free(dfij);
	} // end for (i) for (j)

	return;
}


// Pair Approximation

extern inline double hmca_pa_nn_1x1 (const int *indices, int n_species_0);

extern inline double hmca_pa_nn_2x1 (const int *indices, int n_species_0);

void hmca_pa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;

	int i, j;
	int r0, r1, p0, p1;
	int indices[3];
	double ky, y0, y1;

	memset(dydt, 0, nn_species*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		for (j = 0; j < n_species; ++j)
		{
			ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
			dydt[hmca_sym_id(r0, j, n_species)] -= ky;
			dydt[hmca_sym_id(p0, j, n_species)] += ky;
		}
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];
		dydt[hmca_sym_id(r0, r1, n_species)] -= ky;
		dydt[hmca_sym_id(p0, p1, n_species)] += ky;

		y0 = 0.0;
		y1 = 0.0;
		for (j = 0; j < n_species; ++j)
		{
			y0 += y[hmca_sym_id(r0, j, n_species)];
			y1 += y[hmca_sym_id(r1, j, n_species)];
		}

		for (j = 0; j < n_species; ++j)
		{
			indices[0] = r1;
			indices[1] = r0;
			indices[2] = j;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)] * y[hmca_sym_id(r0, j, n_species)] / (y0 + (y0 == 0.0));
			dydt[hmca_sym_id(r0, j, n_species)] -= ky;
			dydt[hmca_sym_id(p0, j, n_species)] += ky;

			indices[0] = r0;
			indices[1] = r1;
			indices[2] = j;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r1, r0, n_species)] * y[hmca_sym_id(r1, j, n_species)] / (y1 + (y1 == 0.0));
			dydt[hmca_sym_id(r1, j, n_species)] -= ky;
			dydt[hmca_sym_id(p1, j, n_species)] += ky;
		}
	}

	for (i = 0; i < n_species; ++i)
		dydt[hmca_sym_id(i, i, n_species)] += dydt[hmca_sym_id(i, i, n_species)];

	return;
}

void hmca_pa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;

	int i, j, k;
	int r0, r1, p0, p1;
	int indices[3];
	double n, ky, y0, y1;

	memset(dfdy, 0, nn_species*nn_species*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		for (j = 0; j < n_species; ++j)
		{
			// ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
			dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(r0, j, n_species)] -= rates[i];
			dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(r0, j, n_species)] += rates[i];
		}
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		// ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];
		dfdy[nn_species * hmca_sym_id(r0, r1, n_species) + hmca_sym_id(r0, r1, n_species)] -= rates[n_unimol+i];
		dfdy[nn_species * hmca_sym_id(p0, p1, n_species) + hmca_sym_id(r0, r1, n_species)] += rates[n_unimol+i];

		y0 = 0.0;
		y1 = 0.0;
		for (j = 0; j < n_species; ++j)
		{
			y0 += y[hmca_sym_id(r0, j, n_species)];
			y1 += y[hmca_sym_id(r1, j, n_species)];
		}

		for (j = 0; j < n_species; ++j)
		{
			indices[0] = r1;
			indices[1] = r0;
			indices[2] = j;
			n = (*nn) (indices, n_species_0);
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)] * y[hmca_sym_id(r0, j, n_species)] / (y0 + (y0 == 0.0));
			ky = n * rates[n_unimol+i] * (y[hmca_sym_id(r0, j, n_species)] + (y0 == 0.0 && j == r1)) / (y0 + (y0 == 0.0));
			dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(r0, r1, n_species)] -= ky;
			dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(r0, r1, n_species)] += ky;
			ky = n * rates[n_unimol+i] * (y[hmca_sym_id(r0, r1, n_species)] + (y0 == 0.0 && j == r1)) / (y0 + (y0 == 0.0));
			dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(r0, j, n_species)] -= ky;
			dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(r0, j, n_species)] += ky;
			for (k = 0; k < n_species; ++k)
			{
				ky = - n * rates[n_unimol+i] * (y[hmca_sym_id(r0, r1, n_species)] * y[hmca_sym_id(r0, j, n_species)] + (y0 == 0.0 && j == r1 && k == r1)) / (y0 * y0 + (y0 == 0.0));
				dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(r0, k, n_species)] -= ky;
				dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(r0, k, n_species)] += ky;
			}

			indices[0] = r0;
			indices[1] = r1;
			indices[2] = j;
			n = (*nn) (indices, n_species_0);
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r1, r0, n_species)] * y[hmca_sym_id(r1, j, n_species)] / (y1 + (y1 == 0.0));
			ky = n * rates[n_unimol+i] * (y[hmca_sym_id(r1, j, n_species)] + (y1 == 0.0 && j == r0)) / (y1 + (y1 == 0.0));
			dfdy[nn_species * hmca_sym_id(r1, j, n_species) + hmca_sym_id(r1, r0, n_species)] -= ky;
			dfdy[nn_species * hmca_sym_id(p1, j, n_species) + hmca_sym_id(r1, r0, n_species)] += ky;
			ky = n * rates[n_unimol+i] * (y[hmca_sym_id(r1, r0, n_species)] + (y1 == 0.0 && j == r0)) / (y1 + (y1 == 0.0));
			dfdy[nn_species * hmca_sym_id(r1, j, n_species) + hmca_sym_id(r1, j, n_species)] -= ky;
			dfdy[nn_species * hmca_sym_id(p1, j, n_species) + hmca_sym_id(r1, j, n_species)] += ky;
			for (k = 0; k < n_species; ++k)
			{
				ky = - n * rates[n_unimol+i] * (y[hmca_sym_id(r1, r0, n_species)] * y[hmca_sym_id(r1, j, n_species)] + (y1 == 0.0 && j == r0 && k == r0)) / (y1 * y1 + (y1 == 0.0));
				dfdy[nn_species * hmca_sym_id(r1, j, n_species) + hmca_sym_id(r1, k, n_species)] -= ky;
				dfdy[nn_species * hmca_sym_id(p1, j, n_species) + hmca_sym_id(r1, k, n_species)] += ky;
			}
		}
	}

	for (i = 0; i < n_species; ++i) for (j = 0; j < nn_species; ++j)
		dfdy[nn_species * hmca_sym_id(i, i, n_species) + j] += dfdy[nn_species * hmca_sym_id(i, i, n_species) + j];

	return;
}

void hmca_pa_dfdk (
		const double *y,
		double *dfdk,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int n_react = n_unimol+n_bimol;

	int i, j;
	int r0, r1, p0, p1;
	int indices[3];
	double ky, y0, y1;

	memset(dfdk, 0, nn_species*n_react*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		for (j = 0; j < n_species; ++j)
		{
			// ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
			dfdk[n_react * hmca_sym_id(r0, j, n_species) + i] -= y[hmca_sym_id(r0, j, n_species)];
			dfdk[n_react * hmca_sym_id(p0, j, n_species) + i] += y[hmca_sym_id(r0, j, n_species)];
		}
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		// ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];
		dfdk[n_react * hmca_sym_id(r0, r1, n_species) + n_unimol+i] -= y[hmca_sym_id(r0, r1, n_species)];
		dfdk[n_react * hmca_sym_id(p0, p1, n_species) + n_unimol+i] += y[hmca_sym_id(r0, r1, n_species)];

		y0 = 0.0;
		y1 = 0.0;
		for (j = 0; j < n_species; ++j)
		{
			y0 += y[hmca_sym_id(r0, j, n_species)];
			y1 += y[hmca_sym_id(r1, j, n_species)];
		}

		for (j = 0; j < n_species; ++j)
		{
			indices[0] = r1;
			indices[1] = r0;
			indices[2] = j;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)] * y[hmca_sym_id(r0, j, n_species)] / (y0 + (y0 == 0.0));
			ky = (*nn) (indices, n_species_0) * y[hmca_sym_id(r0, r1, n_species)] * y[hmca_sym_id(r0, j, n_species)] / (y0 + (y0 == 0.0));
			dfdk[n_react * hmca_sym_id(r0, j, n_species) + n_unimol+i] -= ky;
			dfdk[n_react * hmca_sym_id(p0, j, n_species) + n_unimol+i] += ky;

			indices[0] = r0;
			indices[1] = r1;
			indices[2] = j;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[hmca_sym_id(r1, r0, n_species)] * y[hmca_sym_id(r1, j, n_species)] / (y1 + (y1 == 0.0));
			ky = (*nn) (indices, n_species_0) * y[hmca_sym_id(r1, r0, n_species)] * y[hmca_sym_id(r1, j, n_species)] / (y1 + (y1 == 0.0));
			dfdk[n_react * hmca_sym_id(r1, j, n_species) + n_unimol+i] -= ky;
			dfdk[n_react * hmca_sym_id(p1, j, n_species) + n_unimol+i] += ky;
		}
	}

	for (i = 0; i < n_species; ++i) for (j = 0; j < n_react; ++j)
		dfdk[n_react * hmca_sym_id(i, i, n_species) + j] += dfdk[n_react * hmca_sym_id(i, i, n_species) + j];

	return;
}


// Half Heterogeneous Pair Approximation

void hmca_hhpa_average (
		const double *y,
		double *_y, double *_ky, double *_yy, double *_kyy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	memset(_y, 0, n_species*n_species*sizeof(double));
	memset(_ky, 0, (n_unimol*n_species+n_bimol)*sizeof(double));
	memset(_yy, 0, n_bimol*n_species*sizeof(double));
	memset(_kyy, 0, n_bimol*n_species*sizeof(double));

	for (int j = 0; j < mesh; ++j)
	{
		const double *kj = rates+n_react*j;
		const double *yj = y+n_species*n_species*j;

		int k, l;
		int r0, r1;
		double y0, y1;

		for (k = 0; k < n_species*n_species; ++k)
		{
			_y[k] += yj[k] * weights[j];
		}

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			for (l = 0; l < n_species; ++l)
				_ky[n_species*k+l] += kj[k] * yj[n_species*r0+l] * weights[j];
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			_ky[n_unimol*n_species+k] += kj[n_unimol+k] * yj[n_species*r0+r1] * weights[j];

			y0 = 0.0;
			y1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				y0 += yj[n_species*r0+l];
				y1 += yj[n_species*r1+l];
			}

			for (l = 0; l < n_species; ++l)
			{
				_yy[n_species*k+l] += yj[n_species*r1+r0] * yj[n_species*r1+l] / (y1 + (y1 == 0.0)) * weights[j];
				_kyy[n_species*k+l] += kj[n_unimol+k] * yj[n_species*r0+r1] * yj[n_species*r0+l] / (y0 + (y0 == 0.0)) * weights[j];
			}
		}
	}

	return;
}

void hmca_hhpa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*n_species*sizeof(double));
	double *_ky = (double*)malloc((n_unimol*n_species+n_bimol)*sizeof(double));
	double *_yy = (double*)malloc(n_bimol*n_species*sizeof(double));
	double *_kyy = (double*)malloc(n_bimol*n_species*sizeof(double));

	hmca_hhpa_average(
			y,
			_y, _ky, _yy, _kyy,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights);

	for (int i = 0; i < mesh; ++i)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*n_species*i;
		double *fi = dydt+n_species*n_species*i;

		int k, l;
		int r0, r1, p0, p1;
		int indices[3];
		double ky, y0, y1;

		memset(fi, 0, n_species*n_species*sizeof(double));

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			p0 = reactions[2*k+1];
			for (l = 0; l < n_species; ++l)
			{
				ky = ki[k] * yi[n_species*r0+l];
				fi[n_species*r0+l] -= ky;
				fi[n_species*p0+l] += ky;
				ky = _ky[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				fi[n_species*l+r0] -= ky;
				fi[n_species*l+p0] += ky;
			}
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			ky = ki[n_unimol+k] * yi[n_species*r0+r1];
			fi[n_species*r0+r1] -= ky;
			fi[n_species*p0+p1] += ky;
			ky = _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
			fi[n_species*r1+r0] -= ky;
			fi[n_species*p1+p0] += ky;

			y0 = 0.0;
			y1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				y0 += yi[n_species*r0+l];
				y1 += yi[n_species*r1+l];
			}

			for (l = 0; l < n_species; ++l)
			{
				indices[0] = r1;
				indices[1] = r0;
				indices[2] = l;
				ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[n_species*r0+r1] * yi[n_species*r0+l] / (y0 + (y0 == 0.0));
				fi[n_species*r0+l] -= ky;
				fi[n_species*p0+l] += ky;
				ky = (*nn) (indices, n_species_0) * _kyy[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				fi[n_species*l+r0] -= ky;
				fi[n_species*l+p0] += ky;

				indices[0] = r0;
				indices[1] = r1;
				indices[2] = l;
				ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l]
					        / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (y1 + (y1 == 0.0));
				fi[n_species*r1+l] -= ky;
				fi[n_species*p1+l] += ky;
				ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1]
					        / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				fi[n_species*l+r1] -= ky;
				fi[n_species*l+p1] += ky;
			}
		}
	}

	free(_y);
	free(_ky);
	free(_yy);
	free(_kyy);

	return;
}

void hmca_hhpa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*n_species*sizeof(double));
	double *_ky = (double*)malloc((n_unimol*n_species+n_bimol)*sizeof(double));
	double *_yy = (double*)malloc(n_bimol*n_species*sizeof(double));
	double *_kyy = (double*)malloc(n_bimol*n_species*sizeof(double));

	hmca_hhpa_average(
			y,
			_y, _ky, _yy, _kyy,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights);

	for (int i = 0; i < mesh; ++i) for (int j = 0; j < mesh; ++j)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*n_species*i;
		const double *kj = rates+n_react*j;
		const double *yj = y+n_species*n_species*j;
		double **dfij = (double**)malloc(n_species*n_species*sizeof(double*));

		int k, l, m;
		int r0, r1, p0, p1;
		int indices[3];
		double n, ky, yi0, yi1, yj0, yj1;

		for (k = 0; k < n_species*n_species; ++k)
		{
			dfij[k] = dfdy + n_species*n_species*mesh*(n_species*n_species*i+k) + n_species*n_species*j;
			memset(dfij[k], 0, n_species*n_species*sizeof(double));
		}

		if (i == j)
		{
			for (k = 0; k < n_unimol; ++k)
			{
				r0 = reactions[2*k+0];
				p0 = reactions[2*k+1];
				for (l = 0; l < n_species; ++l)
				{
					// ky = ki[k] * yi[n_species*r0+l];
					dfij[n_species*r0+l][n_species*r0+l] -= ki[k];
					dfij[n_species*p0+l][n_species*r0+l] += ki[k];
					// ky = _ky[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
					ky = (_ky[n_species*k+l] + (_y[n_species*r0+l] == 0.0 && l == r0) * kj[k]) / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
					dfij[n_species*l+r0][n_species*l+r0] -= ky;
					dfij[n_species*l+p0][n_species*l+r0] += ky;
				}
			}

			for (k = 0; k < n_bimol; ++k)
			{
				r0 = reactions[2*n_unimol+4*k+0];
				r1 = reactions[2*n_unimol+4*k+1];
				p0 = reactions[2*n_unimol+4*k+2];
				p1 = reactions[2*n_unimol+4*k+3];

				// ky = ki[n_unimol+k] * yi[n_species*r0+r1];
				dfij[n_species*r0+r1][n_species*r0+r1] -= ki[n_unimol+k];
				dfij[n_species*p0+p1][n_species*r0+r1] += ki[n_unimol+k];
				// ky = _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
				ky = (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && r0 == r1) * kj[n_unimol+k]) / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
				dfij[n_species*r1+r0][n_species*r1+r0] -= ky;
				dfij[n_species*p1+p0][n_species*r1+r0] += ky;
				
				yi0 = 0.0;
				yi1 = 0.0;
				yj0 = 0.0;
				yj1 = 0.0;
				for (l = 0; l < n_species; ++l)
				{
					yi0 += yi[n_species*r0+l];
					yi1 += yi[n_species*r1+l];
					yj0 += yj[n_species*r0+l];
					yj1 += yj[n_species*r1+l];
				}

				for (l = 0; l < n_species; ++l)
				{
					indices[0] = r1;
					indices[1] = r0;
					indices[2] = l;
					n = (*nn) (indices, n_species_0);
					// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[n_species*r0+r1] * yi[n_species*r0+l] / (yi0 + (yi0 == 0.0));
					ky = n * ki[n_unimol+k] * (yi[n_species*r0+l] + (yi0 == 0.0 && l == r1)) / (yi0 + (yi0 == 0.0));
					dfij[n_species*r0+l][n_species*r0+r1] -= ky;
					dfij[n_species*p0+l][n_species*r0+r1] += ky;
					ky = n * ki[n_unimol+k] * (yi[n_species*r0+r1] + (yi0 == 0.0 && l == r1)) / (yi0 + (yi0 == 0.0));
					dfij[n_species*r0+l][n_species*r0+l] -= ky;
					dfij[n_species*p0+l][n_species*r0+l] += ky;
					for (m = 0; m < n_species; ++m)
					{
						ky = - n * ki[n_unimol+k] * (yi[n_species*r0+r1] * yi[n_species*r0+l] + (yi0 == 0.0 && l == r1 && m == r1)) / (yi0 * yi0 + (yi0 == 0.0));
						dfij[n_species*r0+l][n_species*r0+m] -= ky;
						dfij[n_species*p0+l][n_species*r0+m] += ky;
					}
					// ky = (*nn) (indices, n_species_0) * _kyy[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
					ky = n * (_kyy[n_species*k+l] + (_y[n_species*r0+l] == 0.0 && l == r0) * kj[n_unimol+k] * (yj[n_species*r0+r1] + (yj0 == 0.0 && r0 == r1)) / (yj0 + (yj0 == 0.0)))
						/ (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
					dfij[n_species*l+r0][n_species*l+r0] -= ky;
					dfij[n_species*l+p0][n_species*l+r0] += ky;

					indices[0] = r0;
					indices[1] = r1;
					indices[2] = l;
					n = (*nn) (indices, n_species_0);
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l]
					//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
					ky = n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && r0 == r1) * kj[n_unimol+k])
						* (yi[n_species*r1+l] + (yi1 == 0.0 && l == r0))
						/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
					dfij[n_species*r1+l][n_species*r1+r0] -= ky;
					dfij[n_species*p1+l][n_species*r1+r0] += ky;
					ky = n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && l == r0 && r0 == r1) * kj[n_unimol+k])
						* (yi[n_species*r1+r0] + (yi1 == 0.0 && l == r0))
						/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
					dfij[n_species*r1+l][n_species*r1+l] -= ky;
					dfij[n_species*p1+l][n_species*r1+l] += ky;
					for (m = 0; m < n_species; ++m)
					{
						ky = - n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && m == r0 && r0 == r1) * kj[n_unimol+k])
							* (yi[n_species*r1+r0] * yi[n_species*r1+l] + (yi1 == 0.0 && l == r0 && m == r0))
							/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 * yi1 + (yi1 == 0.0));
						dfij[n_species*r1+l][n_species*r1+m] -= ky;
						dfij[n_species*p1+l][n_species*r1+m] += ky;
					}
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1]
					//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
					ky = n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && l == r0) * kj[n_unimol+k])
						* (_yy[n_species*k+l] + (_y[n_species*r1+l] == 0.0 && l == r1) * (yj[n_species*r1+r0] + (yj1 == 0.0 && r0 == r1)) / (yj1 + (yj1 == 0.0)))
						/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
					dfij[n_species*l+r1][n_species*l+r1] -= ky;
					dfij[n_species*l+p1][n_species*l+r1] += ky;
				}
			}
		} // end if (i == j)

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			p0 = reactions[2*k+1];
			for (l = 0; l < n_species; ++l)
			{
				// ky = ki[k] * yi[n_species*r0+l];
				// ky = _ky[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				ky = kj[k] * weights[j] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0))
					- _ky[n_species*k+l] * yi[n_species*l+r0] * weights[j] / (_y[n_species*r0+l] * _y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				dfij[n_species*l+r0][n_species*r0+l] -= ky;
				dfij[n_species*l+p0][n_species*r0+l] += ky;
			}
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			// ky = ki[n_unimol+k] * yi[n_species*r0+r1];
			// ky = _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
			ky = kj[n_unimol+k] * weights[j] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0))
				- _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * weights[j] / (_y[n_species*r0+r1] * _y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
			dfij[n_species*r1+r0][n_species*r0+r1] -= ky;
			dfij[n_species*p1+p0][n_species*r0+r1] += ky;

			yi0 = 0.0;
			yi1 = 0.0;
			yj0 = 0.0;
			yj1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				yi0 += yi[n_species*r0+l];
				yi1 += yi[n_species*r1+l];
				yj0 += yj[n_species*r0+l];
				yj1 += yj[n_species*r1+l];
			}

			for (l = 0; l < n_species; ++l)
			{
				indices[0] = r1;
				indices[1] = r0;
				indices[2] = l;
				n = (*nn) (indices, n_species_0);
				// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[n_species*r0+r1] * yi[n_species*r0+l] / (yi0 + (yi0 == 0.0));
				// ky = (*nn) (indices, n_species_0) * _kyy[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				ky = n * kj[n_unimol+k] * (yj[n_species*r0+l] + (yj0 == 0.0 && l == r1)) * weights[j] * yi[n_species*l+r0]
					/ (yj0 + (yj0 == 0.0)) / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				dfij[n_species*l+r0][n_species*r0+r1] -= ky;
				dfij[n_species*l+p0][n_species*r0+r1] += ky;
				ky = n * kj[n_unimol+k] * (yj[n_species*r0+r1] + (yj0 == 0.0 && l == r1)) * weights[j] * yi[n_species*l+r0]
					/ (yj0 + (yj0 == 0.0)) / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0))
					- n * _kyy[n_species*k+l] * yi[n_species*l+r0] * weights[j] / (_y[n_species*r0+l] * _y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				dfij[n_species*l+r0][n_species*r0+l] -= ky;
				dfij[n_species*l+p0][n_species*r0+l] += ky;
				for (m = 0; m < n_species; ++m)
				{
					ky = - n * kj[n_unimol+k] * (yj[n_species*r0+r1] * yj[n_species*r0+l] + (yj0 == 0.0 && l == r1 && m == r1)) * weights[j] * yi[n_species*l+r0]
						/ (yj0 * yj0 + (yj0 == 0.0)) / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
					dfij[n_species*l+r0][n_species*r0+m] -= ky;
					dfij[n_species*l+p0][n_species*r0+m] += ky;
				}

				indices[0] = r0;
				indices[1] = r1;
				indices[2] = l;
				n = (*nn) (indices, n_species_0);
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l]
				//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
				ky = n * kj[n_unimol+k] * weights[j] * yi[n_species*r1+r0] * yi[n_species*r1+l]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0))
					- n * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l] * weights[j]
					/ (_y[n_species*r0+r1] * _y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
				dfij[n_species*r1+l][n_species*r0+r1] -= ky;
				dfij[n_species*p1+l][n_species*r0+r1] += ky;
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1]
				//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				ky = n * kj[n_unimol+k] * weights[j]
					* (_yy[n_species*k+l] + (_y[n_species*r0+r1] == 0.0 && r0 == r1) * (yj[n_species*r1+l] + (yj1 == 0.0 && l == r0)) / (yj1 + (yj1 == 0.0))) * yi[n_species*l+r1]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0))
					- n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && r0 == r1) * kj[n_unimol+k])
					* (_yy[n_species*k+l] + (_y[n_species*r0+r1] == 0.0) * (yj[n_species*r1+l] + (yj1 == 0.0 && l == r0)) / (yj1 + (yj1 == 0.0))) * yi[n_species*l+r1] * weights[j]
					/ (_y[n_species*r0+r1] * _y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				dfij[n_species*l+r1][n_species*r0+r1] -= ky;
				dfij[n_species*l+p1][n_species*r0+r1] += ky;
				ky = n * (_ky[n_unimol*n_species+k] + (_y[n_species*r0+r1] == 0.0 && r0 == r1) * kj[n_unimol+k])
					* (yj[n_species*r1+l] + (yj1 == 0.0 && l == r0)) * weights[j] * yi[n_species*l+r1]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yj1 + (yj1 == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				dfij[n_species*l+r1][n_species*r1+r0] -= ky;
				dfij[n_species*l+p1][n_species*r1+r0] += ky;
				ky = n * _ky[n_unimol*n_species+k] * (yj[n_species*r1+r0] + (yj1 == 0.0 && l == r0 && r0 == r1)) * weights[j] * yi[n_species*l+r1]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yj1 + (yj1 == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0))
					- n * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1] * weights[j]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] * _y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				dfij[n_species*l+r1][n_species*r1+l] -= ky;
				dfij[n_species*l+p1][n_species*r1+l] += ky;
				for (m = 0; m < n_species; ++m)
				{
					ky = - n * (_ky[n_unimol*n_species+k] + (_y[n_species*r1+r0] == 0.0 && m == r0 && r0 == r1) * kj[n_unimol+k])
						* (yj[n_species*r1+r0] * yj[n_species*r1+l] + (yj1 == 0.0 && l == r0 && m == r0 && r0 == r1)) * weights[j] * yi[n_species*l+r1]
						/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yj1 * yj1 + (yj1 == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
					dfij[n_species*l+r1][n_species*r1+m] -= ky;
					dfij[n_species*l+p1][n_species*r1+m] += ky;
				}
			}
		}

		free(dfij);
	} // end for (i) for (j)

	free(_y);
	free(_ky);
	free(_yy);
	free(_kyy);

	return;
}

void hmca_hhpa_dfdk (
		const double *y,
		double *dfdk,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(n_species*n_species*sizeof(double));
	double *_ky = (double*)malloc((n_unimol*n_species+n_bimol)*sizeof(double));
	double *_yy = (double*)malloc(n_bimol*n_species*sizeof(double));
	double *_kyy = (double*)malloc(n_bimol*n_species*sizeof(double));

	hmca_hhpa_average(
			y,
			_y, _ky, _yy, _kyy,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights);

	for (int i = 0; i < mesh; ++i) for (int j = 0; j < mesh; ++j)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+n_species*n_species*i;
		const double *kj = rates+n_react*j;
		const double *yj = y+n_species*n_species*j;
		double **dfij = (double**)malloc(n_species*n_species*sizeof(double*));

		int k, l, m;
		int r0, r1, p0, p1;
		int indices[3];
		double ky, yi0, yi1, yj0, yj1;

		for (k = 0; k < n_species*n_species; ++k)
		{
			dfij[k] = dfdk + mesh*n_react*(n_species*n_species*i+k) + n_react*j;
			memset(dfij[k], 0, n_react*sizeof(double));
		}

		if (i == j)
		{
			for (k = 0; k < n_unimol; ++k)
			{
				r0 = reactions[2*k+0];
				p0 = reactions[2*k+1];
				for (l = 0; l < n_species; ++l)
				{
					// ky = ki[k] * yi[n_species*r0+l];
					dfij[n_species*r0+l][k] -= yi[n_species*r0+l];
					dfij[n_species*p0+l][k] += yi[n_species*r0+l];
					// ky = _ky[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				}
			}

			for (k = 0; k < n_bimol; ++k)
			{
				r0 = reactions[2*n_unimol+4*k+0];
				r1 = reactions[2*n_unimol+4*k+1];
				p0 = reactions[2*n_unimol+4*k+2];
				p1 = reactions[2*n_unimol+4*k+3];

				// ky = ki[n_unimol+k] * yi[n_species*r0+r1];
				dfij[n_species*r0+r1][n_unimol+k] -= yi[n_species*r0+r1];
				dfij[n_species*p0+p1][n_unimol+k] += yi[n_species*r0+r1];
				// ky = _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));

				yi0 = 0.0;
				yi1 = 0.0;
				// yj0 = 0.0;
				// yj1 = 0.0;
				for (l = 0; l < n_species; ++l)
				{
					yi0 += yi[n_species*r0+l];
					yi1 += yi[n_species*r1+l];
					// yj0 += yj[n_species*r0+l];
					// yj1 += yj[n_species*r1+l];
				}

				for (l = 0; l < n_species; ++l)
				{
					indices[0] = r1;
					indices[1] = r0;
					indices[2] = l;
					// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[n_species*r0+r1] * yi[n_species*r0+l] / (yi0 + (yi0 == 0.0));
					ky = (*nn) (indices, n_species_0) * yi[n_species*r0+r1] * yi[n_species*r0+l] / (yi0 + (yi0 == 0.0));
					dfij[n_species*r0+l][n_unimol+k] -= ky;
					dfij[n_species*p0+l][n_unimol+k] += ky;
					// ky = (*nn) (indices, n_species_0) * _kyy[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));

					// indices[0] = r0;
					// indices[1] = r1;
					// indices[2] = l;
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l]
					//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1]
					//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				}
			}
		} // end if (i == j)

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			p0 = reactions[2*k+1];
			for (l = 0; l < n_species; ++l)
			{
				// ky = ki[k] * yi[n_species*r0+l];
				// ky = _ky[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				ky = yj[n_species*r0+l] * weights[j] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				dfij[n_species*l+r0][k] -= ky;
				dfij[n_species*l+p0][k] += ky;
			}
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			// ky = ki[n_unimol+k] * yi[n_species*r0+r1];
			// ky = _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
			ky = yj[n_species*r0+r1] * weights[j] * yi[n_species*r1+r0] / (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0));
			dfij[n_species*r1+r0][n_unimol+k] -= ky;
			dfij[n_species*p1+p0][n_unimol+k] += ky;

			yi0 = 0.0;
			yi1 = 0.0;
			yj0 = 0.0;
			yj1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				yi0 += yi[n_species*r0+l];
				yi1 += yi[n_species*r1+l];
				yj0 += yj[n_species*r0+l];
				yj1 += yj[n_species*r1+l];
			}

			for (l = 0; l < n_species; ++l)
			{
				indices[0] = r1;
				indices[1] = r0;
				indices[2] = l;
				// ky = (*nn) (indices, n_species_0) * ki[n_unimol+k] * yi[n_species*r0+r1] * yi[n_species*r0+l] / (yi0 + (yi0 == 0.0));
				// ky = (*nn) (indices, n_species_0) * _kyy[n_species*k+l] * yi[n_species*l+r0] / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				ky = (*nn) (indices, n_species_0) * yj[n_species*r0+r1] * yj[n_species*r0+l] * weights[j] * yi[n_species*l+r0]
					/ (yj0 + (yj0 == 0.0)) / (_y[n_species*r0+l] + (_y[n_species*r0+l] == 0.0));
				dfij[n_species*l+r0][n_unimol+k] -= ky;
				dfij[n_species*l+p0][n_unimol+k] += ky;

				indices[0] = r0;
				indices[1] = r1;
				indices[2] = l;
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * yi[n_species*r1+r0] * yi[n_species*r1+l]
				//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
				ky = (*nn) (indices, n_species_0) * yj[n_species*r0+r1] * weights[j] * yi[n_species*r1+r0] * yi[n_species*r1+l]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (yi1 + (yi1 == 0.0));
				dfij[n_species*r1+l][n_unimol+k] -= ky;
				dfij[n_species*p1+l][n_unimol+k] += ky;
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol*n_species+k] * _yy[n_species*k+l] * yi[n_species*l+r1]
				//	/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				ky = (*nn) (indices, n_species_0) * yj[n_species*r0+r1] * weights[j] * _yy[n_species*k+l] * yi[n_species*l+r1]
					/ (_y[n_species*r0+r1] + (_y[n_species*r0+r1] == 0.0)) / (_y[n_species*r1+l] + (_y[n_species*r1+l] == 0.0));
				dfij[n_species*l+r1][n_unimol+k] -= ky;
				dfij[n_species*l+p1][n_unimol+k] += ky;
			}
		}

		free(dfij);
	} // end for (i) for (j)

	free(_y);
	free(_ky);
	free(_yy);
	free(_kyy);

	return;
}	


// Symmetric Heterogeneous Pair Approximation

void hmca_shpa_average (
		const double *y,
		double *_y, double *_ky,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int n_react = n_unimol+n_bimol;

	memset(_y, 0, nn_species*sizeof(double));
	memset(_ky, 0, n_react*sizeof(double));

	for (int j = 0; j < mesh; ++j)
	{
		const double *kj = rates+n_react*j;
		const double *yj = y+nn_species*j;

		int k, l;
		int r0, r1;
		double yj0;

		for (k = 0; k < nn_species; ++k)
		{
			_y[k] += yj[k] * weights[j];
		}

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];

			yj0 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				yj0 += yj[hmca_sym_id(r0, l, n_species)];
			}

			_ky[k] += kj[k] * yj0 * weights[j];
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			_ky[n_unimol+k] += kj[n_unimol+k] * yj[hmca_sym_id(r0, r1, n_species)] * weights[j];
		}
	}

	return;
}

void hmca_shpa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(nn_species*sizeof(double));
	double *_ky = (double*)malloc(n_react*sizeof(double));

	hmca_shpa_average(
			y,
			_y, _ky,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights
			);

	for (int i = 0; i < mesh; ++i)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+nn_species*i;
		double *fi = dydt+nn_species*i;

		int k, l, m;
		int r0, r1, p0, p1;
		int indices[3];
		double ny, ky, _y0, _y1;

		memset(fi, 0, nn_species*sizeof(double));

		for (k = 0; k < n_unimol; ++k)
		{
			r0 = reactions[2*k+0];
			p0 = reactions[2*k+1];

			_y0 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				_y0 += _y[hmca_sym_id(r0, l, n_species)];
			}

			for (l = 0; l < n_species; ++l)
			{
				ky = ki[k] * yi[hmca_sym_id(r0, l, n_species)];
				fi[hmca_sym_id(r0, l, n_species)] -= ky;
				fi[hmca_sym_id(p0, l, n_species)] += ky;

			//	ny = 0.0;
			//	indices[1] = r0;
			//	indices[2] = l;
			//	for (m = 0; m < n_species; ++m)
			//	{
			//		indices[0] = m;
			//		ny += (*nn) (indices, n_species_0) * _y[hmca_sym_id(r0, m, n_species)];
			//	}

			//	ky = ny * _ky[k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 * _y0 + (_y0 == 0.0));
			//	fi[hmca_sym_id(r0, l, n_species)] -= ky;
			//	fi[hmca_sym_id(p0, l, n_species)] += ky;
			}
		}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			ky = ki[n_unimol+k] * yi[hmca_sym_id(r0, r1, n_species)];
			fi[hmca_sym_id(r0, r1, n_species)] -= ky;
			fi[hmca_sym_id(p0, p1, n_species)] += ky;

			_y0 = 0.0;
			_y1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				_y0 += _y[hmca_sym_id(r0, l, n_species)];
				_y1 += _y[hmca_sym_id(r1, l, n_species)];
			}

			for (l = 0; l < n_species; ++l)
			{
				indices[0] = r1;
				indices[1] = r0;
				indices[2] = l;
				ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 + (_y0 == 0.0));
				fi[hmca_sym_id(r0, l, n_species)] -= ky;
				fi[hmca_sym_id(p0, l, n_species)] += ky;

				indices[0] = r0;
				indices[1] = r1;
				indices[2] = l;
				ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r1, l, n_species)] / (_y1 + (_y1 == 0.0));
				fi[hmca_sym_id(r1, l, n_species)] -= ky;
				fi[hmca_sym_id(p1, l, n_species)] += ky;
			}
		}

		for (k = 0; k < n_species; ++k)
			fi[hmca_sym_id(k, k, n_species)] += fi[hmca_sym_id(k, k, n_species)];
	}

	free(_y);
	free(_ky);

	return;
}

void hmca_shpa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int n_react = n_unimol+n_bimol;

	double *_y = (double*)malloc(nn_species*sizeof(double));
	double *_ky = (double*)malloc(n_react*sizeof(double));

	hmca_shpa_average(
			y,
			_y, _ky,
			n_species_0, n_species_1, n_unimol, n_bimol, mesh,
			reactions, rates, weights
			);

	for (int i = 0; i < mesh; ++i) for (int j = 0; j < mesh; ++j)
	{
		const double *ki = rates+n_react*i;
		const double *yi = y+nn_species*i;
		const double *kj = rates+n_react*j;
		const double *yj = y+nn_species*j;
		double **dfij = (double**)malloc(nn_species*sizeof(double*));

		int k, l, m;
		int r0, r1, p0, p1;
		int indices[3];
		double n, ny, ky, _y0, _y1;

		for (k = 0; k < nn_species; ++k)
		{
			dfij[k] = dfdy + nn_species*mesh*(nn_species*i+k) + nn_species*j;
			memset(dfij[k], 0, nn_species*sizeof(double));
		}

		if (i == j)
		{
			for (k = 0; k < n_unimol; ++k)
			{
				r0 = reactions[2*k+0];
				p0 = reactions[2*k+1];

				_y0 = 0.0;
				for (l = 0; l < n_species; ++l)
				{
					_y0 += _y[hmca_sym_id(r0, l, n_species)];
				}

				for (l = 0; l < n_species; ++l)
				{
					// ky = ki[k] * yi[hmca_sym_id(r0, l, n_species)];
					dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, l, n_species)] -= ki[k];
					dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, l, n_species)] += ki[k];

				//	ny = 0.0;
				//	indices[1] = r0;
				//	indices[2] = l;
				//	for (m = 0; m < n_species; ++m)
				//	{
				//		indices[0] = m;
				//		ny += (*nn) (indices, n_species_0) * _y[hmca_sym_id(r0, m, n_species)];
				//	}

				//	// ky = ny * _ky[k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 * _y0 + (_y0 == 0.0));
				//	indices[0] = l;
				//	ky = (ny * _ky[k] + (_y0 == 0.0) * (*nn) (indices, n_species_0) * kj[k]) / (_y0 * _y0 + (_y0 == 0.0));
				//	dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, l, n_species)] -= ky;
				//	dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, l, n_species)] += ky;
				}
			}

			for (k = 0; k < n_bimol; ++k)
			{
				r0 = reactions[2*n_unimol+4*k+0];
				r1 = reactions[2*n_unimol+4*k+1];
				p0 = reactions[2*n_unimol+4*k+2];
				p1 = reactions[2*n_unimol+4*k+3];

				// ky = ki[n_unimol+k] * yi[hmca_sym_id(r0, r1, n_species)];
				dfij[hmca_sym_id(r0, r1, n_species)][hmca_sym_id(r0, r1, n_species)] -= ki[n_unimol+k];
				dfij[hmca_sym_id(p0, p1, n_species)][hmca_sym_id(r0, r1, n_species)] += ki[n_unimol+k];

				_y0 = 0.0;
				_y1 = 0.0;
				for (l = 0; l < n_species; ++l)
				{
					_y0 += _y[hmca_sym_id(r0, l, n_species)];
					_y1 += _y[hmca_sym_id(r1, l, n_species)];
				}

				for (l = 0; l < n_species; ++l)
				{
					indices[0] = r1;
					indices[1] = r0;
					indices[2] = l;
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 + (_y0 == 0.0));
					ky = (*nn) (indices, n_species_0) * (_ky[n_unimol+k] + (_y0 == 0.0 && l == r1) * kj[n_unimol+k]) / (_y0 + (_y0 == 0.0));
					dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, l, n_species)] -= ky;
					dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, l, n_species)] += ky;

					indices[0] = r0;
					indices[1] = r1;
					indices[2] = l;
					// ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r1, l, n_species)] / (_y1 + (_y1 == 0.0));
					ky = (*nn) (indices, n_species_0) * (_ky[n_unimol+k] + (_y1 == 0.0 && l == r0) * kj[n_unimol+k]) / (_y1 + (_y1 == 0.0));
					dfij[hmca_sym_id(r1, l, n_species)][hmca_sym_id(r1, l, n_species)] -= ky;
					dfij[hmca_sym_id(p1, l, n_species)][hmca_sym_id(r1, l, n_species)] += ky;
				}
			}
		} // end if (i == j)

	//	for (k = 0; k < n_unimol; ++k)
	//	{
	//		r0 = reactions[2*k+0];
	//		p0 = reactions[2*k+1];

	//		_y0 = 0.0;
	//		for (l = 0; l < n_species; ++l)
	//		{
	//			_y0 += _y[hmca_sym_id(r0, l, n_species)];
	//		}

	//		for (l = 0; l < n_species; ++l)
	//		{
	//			// ky = ki[k] * yi[hmca_sym_id(r0, l, n_species)];

	//			ny = 0.0;
	//			indices[1] = r0;
	//			indices[2] = l;
	//			for (m = 0; m < n_species; ++m)
	//			{
	//				indices[0] = m;
	//				ny += (*nn) (indices, n_species_0) * _y[hmca_sym_id(r0, m, n_species)];
	//			}

	//			// ky = ny * _ky[k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 * _y0 + (_y0 == 0.0));
	//			ky = ny * kj[k] * weights[j] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 * _y0 + (_y0 == 0.0))
	//				- 2.0 * ny * _ky[k] * yi[hmca_sym_id(r0, l, n_species)] * weights[j] / (_y0 * _y0 * _y0 + (_y0 == 0.0));
	//			for (m = 0; m < n_species; ++m)
	//			{
	//				dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, m, n_species)] -= ky;
	//				dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, m, n_species)] += ky;
	//			}
	//			for (m = 0; m < n_species; ++m)
	//			{
	//				indices[0] = m;
	//				ky = (*nn) (indices, n_species_0) * weights[j] * _ky[k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 * _y0 + (_y0 == 0.0));
	//				dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, m, n_species)] -= ky;
	//				dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, m, n_species)] += ky;
	//			}
	//		}
	//	}

		for (k = 0; k < n_bimol; ++k)
		{
			r0 = reactions[2*n_unimol+4*k+0];
			r1 = reactions[2*n_unimol+4*k+1];
			p0 = reactions[2*n_unimol+4*k+2];
			p1 = reactions[2*n_unimol+4*k+3];

			// ky = ki[n_unimol+k] * yi[hmca_sym_id(r0, r1, n_species)];

			_y0 = 0.0;
			_y1 = 0.0;
			for (l = 0; l < n_species; ++l)
			{
				_y0 += _y[hmca_sym_id(r0, l, n_species)];
				_y1 += _y[hmca_sym_id(r1, l, n_species)];
			}

			for (l = 0; l < n_species; ++l)
			{
				indices[0] = r1;
				indices[1] = r0;
				indices[2] = l;
				n = (*nn) (indices, n_species_0);
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 + (_y0 == 0.0));
				ky = n * kj[n_unimol+k] * weights[j] * yi[hmca_sym_id(r0, l, n_species)] / (_y0 + (_y0 == 0.0));
				dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, r1, n_species)] -= ky;
				dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, r1, n_species)] += ky;
				for (m = 0; m < n_species; ++m)
				{
					ky = - n * _ky[n_unimol+k] * yi[hmca_sym_id(r0, l, n_species)] * weights[j] / (_y0 * _y0 + (_y0 == 0.0));
					dfij[hmca_sym_id(r0, l, n_species)][hmca_sym_id(r0, m, n_species)] -= ky;
					dfij[hmca_sym_id(p0, l, n_species)][hmca_sym_id(r0, m, n_species)] += ky;
				}

				indices[0] = r0;
				indices[1] = r1;
				indices[2] = l;
				n = (*nn) (indices, n_species_0);
				// ky = (*nn) (indices, n_species_0) * _ky[n_unimol+k] * yi[hmca_sym_id(r1, l, n_species)] / (_y1 + (_y1 == 0.0));
				ky = n * kj[n_unimol+k] * weights[j] * yi[hmca_sym_id(r1, l, n_species)] / (_y1 + (_y1 == 0.0));
				dfij[hmca_sym_id(r1, l, n_species)][hmca_sym_id(r1, r0, n_species)] -= ky;
				dfij[hmca_sym_id(p1, l, n_species)][hmca_sym_id(r1, r0, n_species)] += ky;
				for (m = 0; m < n_species; ++m)
				{
					ky = - n * _ky[n_unimol+k] * yi[hmca_sym_id(r1, l, n_species)] * weights[j] / (_y1 * _y1 + (_y1 == 0.0));
					dfij[hmca_sym_id(r1, l, n_species)][hmca_sym_id(r1, m, n_species)] -= ky;
					dfij[hmca_sym_id(p1, l, n_species)][hmca_sym_id(r1, m, n_species)] += ky;
				}
			}
		}

		for (k = 0; k < n_species; ++k) for (l = 0; l < nn_species; ++l)
			dfij[hmca_sym_id(k, k, n_species)][l] += dfij[hmca_sym_id(k, k, n_species)][l];

		free(dfij);
	} // end for (i) for (j)

	free(_y);
	free(_ky);

	return;
}

void hmca_shpa_dfdk (
		const double *y,
		double *dfdk,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
		const int *reactions, const double *rates, const double *weights, hmca_nn nn
		)
{
	return;
}


// Machine Learning Moment Closure

void hmca_mlmc_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		hmca_mc closure, hmca_mc deriv, void *params
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;

	int i, j;
	int r0, r1, p0, p1;
	int indices[3];
	double ky;

	memset(dydt, 0, nn_species*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		for (j = 0; j < n_species; ++j)
		{
			ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
			dydt[hmca_sym_id(r0, j, n_species)] -= ky;
			dydt[hmca_sym_id(p0, j, n_species)] += ky;
		}
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];
		dydt[hmca_sym_id(r0, r1, n_species)] -= ky;
		dydt[hmca_sym_id(p0, p1, n_species)] += ky;

		for (j = 0; j < n_species; ++j)
		{
			indices[0] = r1;
			indices[1] = r0;
			indices[2] = j;
			ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			dydt[hmca_sym_id(r0, j, n_species)] -= ky;
			dydt[hmca_sym_id(p0, j, n_species)] += ky;

			indices[0] = r0;
			indices[1] = r1;
			indices[2] = j;
			ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			dydt[hmca_sym_id(r1, j, n_species)] -= ky;
			dydt[hmca_sym_id(p1, j, n_species)] += ky;
		}
	}

	for (i = 0; i < n_species; ++i)
		dydt[hmca_sym_id(i, i, n_species)] += dydt[hmca_sym_id(i, i, n_species)];

	return;
}

void hmca_mlmc_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		hmca_mc closure, hmca_mc deriv, void *params
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;

	int i, j, k, l;
	int r0, r1, p0, p1;
	int indices[5];
	double ky;

	memset(dfdy, 0, nn_species*nn_species*sizeof(double));

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		for (j = 0; j < n_species; ++j)
		{
			// ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
			dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(r0, j, n_species)] -= rates[i];
			dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(r0, j, n_species)] += rates[i];
		}
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		// ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];
		dfdy[nn_species * hmca_sym_id(r0, r1, n_species) + hmca_sym_id(r0, r1, n_species)] -= rates[n_unimol+i];
		dfdy[nn_species * hmca_sym_id(p0, p1, n_species) + hmca_sym_id(r0, r1, n_species)] += rates[n_unimol+i];
		
		for (j = 0; j < n_species; ++j)
		{
			indices[0] = r1;
			indices[1] = r0;
			indices[2] = j;
			// ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			for (k = 0; k < n_species; ++k)
			{
				indices[3] = k;
				for (l = k; l < n_species; ++l)
				{
					indices[4] = l;
					ky = rates[n_unimol+i] * (*deriv) (indices, y, n_species_0, n_species_1, params);
					dfdy[nn_species * hmca_sym_id(r0, j, n_species) + hmca_sym_id(k, l, n_species)] -= ky;
					dfdy[nn_species * hmca_sym_id(p0, j, n_species) + hmca_sym_id(k, l, n_species)] += ky;
				}
			}

			indices[0] = r0;
			indices[1] = r1;
			indices[2] = j;
			// ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			for (k = 0; k < n_species; ++k)
			{
				indices[3] = k;
				for (l = k; l < n_species; ++l)
				{
					indices[4] = l;
					ky = rates[n_unimol+i] * (*deriv) (indices, y, n_species_0, n_species_1, params);
					dfdy[nn_species * hmca_sym_id(r1, j, n_species) + hmca_sym_id(k, l, n_species)] -= ky;
					dfdy[nn_species * hmca_sym_id(p1, j, n_species) + hmca_sym_id(k, l, n_species)] += ky;
				}
			}
		}
	}

	for (i = 0; i < n_species; ++i) for (j = 0; j < nn_species; ++j)
		dfdy[nn_species * hmca_sym_id(i, i, n_species) + j] += dfdy[nn_species * hmca_sym_id(i, i, n_species) + j];

	return;
}

void hmca_mlmc_dfdz (
		const double *y,
		double *dfdz,
		int n_species_0, int n_species_1, int n_unimol, int n_bimol,
		const int *reactions, const double *rates,
		hmca_mc closure, hmca_mc deriv, void *params
		)
{
	int n_species = n_species_0+n_species_1;
	int nn_species = n_species*(n_species+1)/2;
	int nnn_species = n_species*nn_species;

	int i, j, k, l;
	int r0, r1, p0, p1;
	int index;
	double ky;

	memset(dfdz, 0, nn_species*nnn_species*sizeof(double));

	// for (i = 0; i < n_unimol; ++i)
	// {
	//	r0 = reactions[2*i+0];
	//	p0 = reactions[2*i+1];
	//	for (j = 0; j < n_species; ++j)
	//	{
	//		ky = rates[i] * y[hmca_sym_id(r0, j, n_species)];
	//	}
	// }

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		// ky = rates[n_unimol+i] * y[hmca_sym_id(r0, r1, n_species)];

		for (j = 0; j < n_species; ++j)
		{
			// indices[0] = r1;
			// indices[1] = r0;
			// indices[2] = j;
			index = nn_species*r0+hmca_sym_id(r1, j, n_species);
			// ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			dfdz[nnn_species * hmca_sym_id(r0, j, n_species) + index] -= rates[n_unimol+i];
			dfdz[nnn_species * hmca_sym_id(p0, j, n_species) + index] += rates[n_unimol+i];

			// indices[0] = r0;
			// indices[1] = r1;
			// indices[2] = j;
			index = nn_species*r1+hmca_sym_id(r0, j, n_species);
			// ky = rates[n_unimol+i] * (*closure) (indices, y, n_species_0, n_species_1, params);
			dfdz[nnn_species * hmca_sym_id(r1, j, n_species) + index] -= rates[n_unimol+i];
			dfdz[nnn_species * hmca_sym_id(p1, j, n_species) + index] += rates[n_unimol+i];
		}
	}

	for (i = 0; i < n_species; ++i) for (j = 0; j < nnn_species; ++j)
		dfdz[nnn_species * hmca_sym_id(i, i, n_species) + j] += dfdz[nnn_species * hmca_sym_id(i, i, n_species) + j];

	return;
}


// Select Pair Approximation

extern inline double hmca_spa_nn_1x1 (const int *indices, int n_species_0);

extern inline double hmca_spa_nn_2x1 (const int *indices, int n_species_0);

void hmca_spa_func (
		const double *y,
		double *dydt,
		int n_species_0, int n_species_1, int n_pairs, int n_unimol, int n_bimol,
		const int *pairs, const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_dof = n_species+n_pairs;

	int i, j, k, l, p, q;
	int r0, r1, p0, p1, s0, s1;
	int indices[3];
	double ky;

	memset(dydt, 0, n_dof*sizeof(double));

	indices[2] = -1;

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		ky = rates[i] * y[r0];
		dydt[r0] -= ky;
		dydt[p0] += ky;
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		for (j = 0; j < n_pairs; ++j)
		{
			if ((r0 == pairs[2*j+0] && r1 == pairs[2*j+1]) || (r0 == pairs[2*j+1] && r1 == pairs[2*j+0]))
				break;
		}
		if (j < n_pairs)
		{
			indices[0] = r0;
			indices[1] = r1;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[n_species+j];
			dydt[r0] -= ky;
			dydt[p0] += ky;

			indices[0] = r1;
			indices[1] = r0;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[n_species+j];
			dydt[r1] -= ky;
			dydt[p1] += ky;
		}
		else
		{
			indices[0] = r0;
			indices[1] = r1;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
			dydt[r0] -= ky;
			dydt[p0] += ky;

			indices[0] = r1;
			indices[1] = r0;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
			dydt[r1] -= ky;
			dydt[p1] += ky;
		}
	} // end for (i)

	for (i = 0; i < n_pairs; ++i)
	{
		s0 = pairs[2*i+0];
		s1 = pairs[2*i+1];

		for (j = 0; j < n_unimol; ++j)
		{
			r0 = reactions[2*j+0];
			p0 = reactions[2*j+1];
			dydt[n_species+i] -= ((s0 == r0) + (s1 == r0)) * rates[j] * y[n_species+i];
			for (k = 0; k < 2; ++k) if (pairs[2*i+k] == p0)
			{
				for (l = 0; l < n_pairs; ++l)
				{
					if ((r0 == pairs[2*l+0] && pairs[2*i+(k+1)%2] == pairs[2*l+1]) || (r0 == pairs[2*l+1] && pairs[2*i+(k+1)%2] == pairs[2*l+0]))
						break;
				}
				if (l < n_pairs)
					dydt[n_species+i] += rates[j] * y[n_species+l];
				else
					dydt[n_species+i] += rates[j] * y[r0] * y[pairs[2*i+(k+1)%2]];
			}
		}

		for (j = 0; j < n_bimol; ++j)
		{
			r0 = reactions[2*n_unimol+4*j+0];
			r1 = reactions[2*n_unimol+4*j+1];
			p0 = reactions[2*n_unimol+4*j+2];
			p1 = reactions[2*n_unimol+4*j+3];

			dydt[n_species+i] -= ((s0 == r0 && s1 == r1) + (s0 == r1 && s1 == r0)) * rates[n_unimol+j] * y[n_species+i];
			for (k = 0; k < n_pairs; ++k)
			{
				if ((r0 == pairs[2*k+0] && r1 == pairs[2*k+1]) || (r0 == pairs[2*k+1] && r1 == pairs[2*k+0]))
					break;
			}
			if (k < n_pairs)
				dydt[n_species+i] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[n_species+k];
			else
				dydt[n_species+i] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[r0] * y[r1];

			for (k = 0; k < 2; ++k) for (l = 0; l < 2; ++l)
			{
				if (pairs[2*i+k] == reactions[2*n_unimol+4*j+l])
				{
					indices[0] = reactions[2*n_unimol+4*j+(l+1)%2];
					indices[1] = reactions[2*n_unimol+4*j+l];
					indices[2] = pairs[2*i+(k+1)%2];
					for (p = 0; p < n_pairs; ++p)
					{
						if ((indices[0] == pairs[2*p+0] && indices[1] == pairs[2*p+1]) || (indices[0] == pairs[2*p+1] && indices[1] == pairs[2*p+0]))
							break;
					}
					for (q = 0; q < n_pairs; ++q)
					{
						if ((indices[1] == pairs[2*q+0] && indices[2] == pairs[2*q+1]) || (indices[1] == pairs[2*q+1] && indices[2] == pairs[2*q+0]))
							break;
					}
					if (p < n_pairs)
					{
						if (q < n_pairs)
							dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
						else
							dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[indices[2]];
					}
					else
					{
						if (q < n_pairs)
							dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[n_species+q];
						else
							dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[indices[1]] * y[indices[2]];
					}
				}

				if (pairs[2*i+k] == reactions[2*n_unimol+4*j+2+l])
				{
					indices[0] = reactions[2*n_unimol+4*j+(l+1)%2];
					indices[1] = reactions[2*n_unimol+4*j+l];
					indices[2] = pairs[2*i+(k+1)%2];
					for (p = 0; p < n_pairs; ++p)
					{
						if ((indices[0] == pairs[2*p+0] && indices[1] == pairs[2*p+1]) || (indices[0] == pairs[2*p+1] && indices[1] == pairs[2*p+0]))
							break;
					}
					for (q = 0; q < n_pairs; ++q)
					{
						if ((indices[1] == pairs[2*q+0] && indices[2] == pairs[2*q+1]) || (indices[1] == pairs[2*q+1] && indices[2] == pairs[2*q+0]))
							break;
					}
					if (p < n_pairs)
					{
						if (q < n_pairs)
							dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
						else
							dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[indices[2]];
					}
					else
					{
						if (q < n_pairs)
							dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[n_species+q];
						else
							dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[indices[1]] * y[indices[2]];
					}
				}
			}
		} // end for (j)
	} // end for (i)

	return;
}

void hmca_spa_jac (
		const double *y,
		double *dfdy,
		int n_species_0, int n_species_1, int n_pairs, int n_unimol, int n_bimol,
		const int *pairs, const int *reactions, const double *rates, hmca_nn nn
		)
{
	int n_species = n_species_0+n_species_1;
	int n_dof = n_species+n_pairs;

	int i, j, k, l, p, q;
	int r0, r1, p0, p1, s0, s1;
	int indices[3];
	double n, ky;

	memset(dfdy, 0, n_dof*n_dof*sizeof(double));

	indices[2] = -1;

	for (i = 0; i < n_unimol; ++i)
	{
		r0 = reactions[2*i+0];
		p0 = reactions[2*i+1];
		// ky = rates[i] * y[r0];
		dfdy[n_dof*r0+r0] -= rates[i];
		dfdy[n_dof*p0+r0] += rates[i];
	}

	for (i = 0; i < n_bimol; ++i)
	{
		r0 = reactions[2*n_unimol+4*i+0];
		r1 = reactions[2*n_unimol+4*i+1];
		p0 = reactions[2*n_unimol+4*i+2];
		p1 = reactions[2*n_unimol+4*i+3];

		for (j = 0; j < n_pairs; ++j)
		{
			if ((r0 == pairs[2*j+0] && r1 == pairs[2*j+1]) || (r0 == pairs[2*j+1] && r1 == pairs[2*j+0]))
				break;
		}
		if (j < n_pairs)
		{
			indices[0] = r0;
			indices[1] = r1;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[n_species+j];
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i];
			dfdy[n_dof*r0+(n_species+j)] -= ky;
			dfdy[n_dof*p0+(n_species+j)] += ky;
			
			indices[0] = r1;
			indices[1] = r0;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[n_species+j];
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i];
			dfdy[n_dof*r1+(n_species+j)] -= ky;
			dfdy[n_dof*p1+(n_species+j)] += ky;
		}
		else
		{
			indices[0] = r0;
			indices[1] = r1;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r1];
			dfdy[n_dof*r0+r0] -= ky;
			dfdy[n_dof*p0+r0] += ky;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0];
			dfdy[n_dof*r0+r1] -= ky;
			dfdy[n_dof*p0+r1] += ky;
			
			indices[0] = r1;
			indices[1] = r0;
			// ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0] * y[r1];
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r1];
			dfdy[n_dof*r1+r0] -= ky;
			dfdy[n_dof*p1+r0] += ky;
			ky = (*nn) (indices, n_species_0) * rates[n_unimol+i] * y[r0];
			dfdy[n_dof*r1+r1] -= ky;
			dfdy[n_dof*p1+r1] += ky;
		}
	}

	for (i = 0; i < n_pairs; ++i)
	{
		s0 = pairs[2*i+0];
		s1 = pairs[2*i+1];

		for (j = 0; j < n_unimol; ++j)
		{
			r0 = reactions[2*j+0];
			p0 = reactions[2*j+1];
			// dydt[n_species+i] -= ((s0 == r0) + (s1 == r0)) * rates[j] * y[n_species+i];
			dfdy[n_dof*(n_species+i)+(n_species+i)] -= ((s0 == r0) + (s1 == r0)) * rates[j];
			for (k = 0; k < 2; ++k) if (pairs[2*i+k] == p0)
			{
				for (l = 0; l < n_pairs; ++l)
				{
					if ((r0 == pairs[2*l+0] && pairs[2*i+(k+1)%2] == pairs[2*l+1]) || (r0 == pairs[2*l+1] && pairs[2*i+(k+1)%2] == pairs[2*l+0]))
						break;
				}
				if (l < n_pairs)
				{
					// dydt[n_species+i] += rates[j] * y[n_species+l];
					dfdy[n_dof*(n_species+i)+(n_species+l)] += rates[j];
				}
				else
				{
					// dydt[n_species+i] += rates[j] * y[r0] * y[pairs[2*i+(k+1)%2]];
					dfdy[n_dof*(n_species+i)+r0] += rates[j] * y[pairs[2*i+(k+1)%2]];
					dfdy[n_dof*(n_species+i)+pairs[2*i+(k+1)%2]] += rates[j] * y[r0];
				}
			}
		}

		for (j = 0; j < n_bimol; ++j)
		{
			r0 = reactions[2*n_unimol+4*j+0];
			r1 = reactions[2*n_unimol+4*j+1];
			p0 = reactions[2*n_unimol+4*j+2];
			p1 = reactions[2*n_unimol+4*j+3];

			// dydt[n_species+i] -= ((s0 == r0 && s1 == r1) + (s0 == r1 && s1 == r0)) * rates[n_unimol+j] * y[n_species+i];
			dfdy[n_dof*(n_species+i)+(n_species+i)] -= ((s0 == r0 && s1 == r1) + (s0 == r1 && s1 == r0)) * rates[n_unimol+j];
			for (k = 0; k < n_pairs; ++k)
			{
				if ((r0 == pairs[2*k+0] && r1 == pairs[2*k+1]) || (r0 == pairs[2*k+1] && r1 == pairs[2*k+0]))
					break;
			}
			if (k < n_pairs)
			{
				// dydt[n_species+i] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[n_species+k];
				dfdy[n_dof*(n_species+i)+(n_species+k)] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j];
			}
			else
			{
				// dydt[n_species+i] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[r0] * y[r1];
				dfdy[n_dof*(n_species+i)+r0] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[r1];
				dfdy[n_dof*(n_species+i)+r1] += ((s0 == p0 && s1 == p1) + (s0 == p1 && s1 == p0)) * rates[n_unimol+j] * y[r0];
			}

			for (k = 0; k < 2; ++k) for (l = 0; l < 2; ++l)
			{
				if (pairs[2*i+k] == reactions[2*n_unimol+4*j+l])
				{
					indices[0] = reactions[2*n_unimol+4*j+(l+1)%2];
					indices[1] = reactions[2*n_unimol+4*j+l];
					indices[2] = pairs[2*i+(k+1)%2];
					n = (*nn) (indices, n_species_0);
					for (p = 0; p < n_pairs; ++p)
					{
						if ((indices[0] == pairs[2*p+0] && indices[1] == pairs[2*p+1]) || (indices[0] == pairs[2*p+1] && indices[1] == pairs[2*p+0]))
							break;
					}
					for (q = 0; q < n_pairs; ++q)
					{
						if ((indices[1] == pairs[2*q+0] && indices[2] == pairs[2*q+1]) || (indices[1] == pairs[2*q+1] && indices[2] == pairs[2*q+0]))
							break;
					}
					if (p < n_pairs)
					{
						if (q < n_pairs)
						{
							// dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+(n_species+p)] -= n * rates[n_unimol+j] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+(n_species+q)] -= n * rates[n_unimol+j] * y[n_species+p] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+indices[1]] -= - n * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] * y[indices[1]] + (y[indices[1]] == 0.0));
						}
						else
						{
							// dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+(n_species+p)] -= n * rates[n_unimol+j] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[2]] -= n * rates[n_unimol+j] * y[n_species+p];
						}
					}
					else
					{
						if (q < n_pairs)
						{
							// dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[n_species+q];
							dfdy[n_dof*(n_species+i)+indices[0]] -= n * rates[n_unimol+j] * y[n_species+q];
							dfdy[n_dof*(n_species+i)+(n_species+q)] -= n * rates[n_unimol+j] * y[indices[0]];
						}
						else
						{
							// dydt[n_species+i] -= (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[indices[1]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[0]] -= n * rates[n_unimol+j] * y[indices[1]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[1]] -= n * rates[n_unimol+j] * y[indices[0]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[2]] -= n * rates[n_unimol+j] * y[indices[0]] * y[indices[1]];
						}
					}
				}

				if (pairs[2*i+k] == reactions[2*n_unimol+4*j+2+l])
				{
					indices[0] = reactions[2*n_unimol+4*j+(l+1)%2];
					indices[1] = reactions[2*n_unimol+4*j+l];
					indices[2] = pairs[2*i+(k+1)%2];
					n = (*nn) (indices, n_species_0);
					for (p = 0; p < n_pairs; ++p)
					{
						if ((indices[0] == pairs[2*p+0] && indices[1] == pairs[2*p+1]) || (indices[0] == pairs[2*p+1] && indices[1] == pairs[2*p+0]))
							break;
					}
					for (q = 0; q < n_pairs; ++q)
					{
						if ((indices[1] == pairs[2*q+0] && indices[2] == pairs[2*q+1]) || (indices[1] == pairs[2*q+1] && indices[2] == pairs[2*q+0]))
							break;
					}
					if (p < n_pairs)
					{
						if (q < n_pairs)
						{
							// dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+(n_species+p)] += n * rates[n_unimol+j] * y[n_species+q] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+(n_species+q)] += n * rates[n_unimol+j] * y[n_species+p] / (y[indices[1]] + (y[indices[1]] == 0.0));
							dfdy[n_dof*(n_species+i)+indices[1]] += - n * rates[n_unimol+j] * y[n_species+p] * y[n_species+q] / (y[indices[1]] * y[indices[1]] + (y[indices[1]] == 0.0));
						}
						else
						{
							// dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[n_species+p] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+(n_species+p)] += n * rates[n_unimol+j] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[2]] += n * rates[n_unimol+j] * y[n_species+p];
						}
					}
					else
					{
						if (q < n_pairs)
						{
							// dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[n_species+q];
							dfdy[n_dof*(n_species+i)+indices[0]] += n * rates[n_unimol+j] * y[n_species+q];
							dfdy[n_dof*(n_species+i)+(n_species+q)] += n * rates[n_unimol+j] * y[indices[0]];
						}
						else
						{
							// dydt[n_species+i] += (*nn) (indices, n_species_0) * rates[n_unimol+j] * y[indices[0]] * y[indices[1]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[0]] += n * rates[n_unimol+j] * y[indices[1]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[1]] += n * rates[n_unimol+j] * y[indices[0]] * y[indices[2]];
							dfdy[n_dof*(n_species+i)+indices[2]] += n * rates[n_unimol+j] * y[indices[0]] * y[indices[1]];
						}
					}
				}
			} // end for (k) for (l)
		} // end for (j)
	} // end for (i)

	return;
}

