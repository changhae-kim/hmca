# Heterogeneous Moment Closure Approximation (HMCA)

## Overview

This code defines a number of functions that give the right-hand sides of the kinetic equations in various methods: MF, PA, HMF, HHPA, SHPA, SPA, and MLMC.

### Example

To compile `test_pa.c`, run

    gcc test_pa.c hmca.c -lm

## Common Parameters

There are some parameters that the functions require in common.

    int n_species_0 = number of species on the first type of sites
    int n_species_1 = number of species on the second type of sites, 0 on a 1x1 lattice
    int n_unimol    = number of unimolecular reactions
    int n_bimol     = number of bimolecular reactions
    int *reactions  = array of the 2*n_unimol+4*n_bimol reactants and products
    double *rates   = array of the n_unimol+n_bimol rate constants
    hmca_nn nn      = function returning the number of nearest neighbors

Two parameters are needed to give the number of species, because 2x1 or 2x2 lattices have two types of sites.

Assign each of the `n_species_0` species to an integer index between `0` and `n_species_0-1`,
and each of the `n_species_1` species to an integer index between `n_species_0` and `n_species_0+n_species_1-1`.

List the unimolecular reactions and then the bimolecular reactions - i.e. reactant, product, reactant, product, . . . reactant 1, reactant 2, product 1, product 2, . . .

The number of nearest neighbors must be given as a function, because it can depend on the spatial arrangement of the sites in question.

    double (*hmca_nn) (const int *indices, int n_species_0)

There are pre-defined functions.

    double hmca_mf_nn_1x1 (const int *indices, int n_species_0)
    double hmca_mf_nn_2x1 (const int *indices, int n_species_0)
    double hmca_pa_nn_1x1 (const int *indices, int n_species_0)
    double hmca_pa_nn_2x1 (const int *indices, int n_species_0)
    double hmca_spa_nn_1x1 (const int *indices, int n_species_0)
    double hmca_spa_nn_2x1 (const int *indices, int n_species_0)

### Example

Consider the Langmuir-Hinshelwood mechanism:

        O --> A        k1 = 1.0
        O --> B        k2 = 1.0
    A + B --> O + O    k3 = 10.0

The code would look like:

    enum {O, A, B};
    int n_species_0 = 3;
    int n_species_1 = 0;
    int n_unimol    = 2;
    int n_bimol     = 1;
    int reactions[] = {O,A, O,B, A,B,O,O};
    double rates[]  = {1.0, 1.0, 10.0};
    
    double *y       = (double*)malloc((n_species_0+n_species_1)*sizeof(double));
    double *dydt    = (double*)malloc((n_species_0+n_species_1)*sizeof(double));
    
    hmca_mf_func(
        y,
        dydt,
        n_species_0, n_species_1, n_unimol, n_bimol,
        reactions, rates, hmca_mf_nn_1x1
        );

## Common Parameters in Heterogeneous Methods

There are some more parameters that the heterogeneous methods require in common,
and one that has a different definition.

    int mesh        = number of points to sample in the rate constant space
    double *rates   = array of the mesh*(n_unimol+n_bimol) rate constants
    double *weights = array of the mesh weights

There are pre-defined functions to set `rates` and `weights` using typical distributions:

    void hmca_lognorm_set ( . . . );
    void hmca_logexp_set ( . . . );
    void hmca_logsech_set ( . . . );
    void hmca_logpoisson2_set (
        const double *logk0, const double *dlogk,
        double *rates, double *weights,
        int n_unimol, int n_bimol, int mesh, double bound
        );

Most of the parameters are as described above.

    const double *logk0 = input,  array of the "zero-point" rate constants in log space
    const double *dlogk = input,  array of the rate constant spreads in log space - i.e. sigma in log-normal distribution
    double *rates       = output, array of the mesh*(n_unimol+n_bimol) rate constants
    double *weights     = output, array of the mesh weights
    double bound        = param,  number of standard deviations to scan

The recommended combinations of `mesh` and `bound` depend on the distribution.

For the log-normal distribution, we recommend:

    double bound = sigma+3.0;
    int mesh = (int)(3.0*(sigma+3.0)+1.0);

For the log-exponential distribution, we recommend:

    double bound = 6.0;
    int mesh = (int)(6.0*(lambda+1.0)+1.0);

For the log-Poisson distribution (n = 2), we recommend:

    double bound = 10.0;
    int mesh = (int)(10.0*(lambda+1.0)+1.0);

### Example

Suppose that the Langmuir-Hinshelwood mechanism has a static disorder in the last step:

        O --> A        mu = 0.0
        O --> B        mu = 0.0
    A + B --> O + O    mu = 1.0    sigma = 5.0

The code would look like:

    int mesh        = (int)(3.0*(5.0+3.0)+1.0);
    double bound    = 5.0+3.0;
    double logk0[]  = {0.0, 0.0, 1.0};
    double dlogk[]  = {0.0, 0.0, 5.0};
    double *rates   = (double*)malloc(mesh*(n_unimol+n_bimol)*sizeof(double));
    double *weights = (double*)malloc(mesh*sizeof(double));

    hmca_lognorm_set (
        logk0, dlogk,
        rates, weights,
        n_unimol, n_bimol, mesh, bound
        );

## Mean-Field Approximation

    void hmca_mf_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol,
        const int *reactions, const double *rates, hmca_nn nn
        )

Most of the parameters are as described above.

    double *y    = input,  coverages of n_species_0+n_species_1 species
    double *dydt = output, rates

## Heterogeneous Mean-Field Approximation

    void hmca_hmf_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
        const int *reactions, const double *rates, const double *weights, hmca_nn nn
        );

Most of the parameters are as described above.

    double *y    = input,  coverages of mesh*(n_species_0+n_species_1) species
    double *dydt = output, rates

## Pair Approximation

    void hmca_pa_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol,
        const int *reactions, const double *rates, hmca_nn nn
        )

Most of the parameters are as described above.

    double *y    = input,  coverages of n_species*(n_species+1)/2 pairs, n_species = n_species_0+n_species_1
    double *dydt = output, rates
    
Due to the symmetry, there are only `n_species*(n_species+1)/2` distinct pairs.

## Half Heterogeneous Pair Approximation

    void hmca_hhpa_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
        const int *reactions, const double *rates, const double *weights, hmca_nn nn
        )

Most of the parameters are as described above.

    double *y    = input,  coverages of mesh*n_species*n_species pairs, n_species = n_species_0+n_species_1
    double *dydt = output, rates

Due to the broken symmetry, there are `n_species*n_species` distinct pairs per each of the `mesh` points.
