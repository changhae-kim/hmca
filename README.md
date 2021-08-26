# Heterogeneous Moment Closure Approximation (HMCA)

## Overview

This code defines a number of functions that give the right-hand sides of the kinetic equations in various methods: MF, PA, HMF, HHPA, SHPA, SPA, and MLMC.

## Common Parameters

There are some parameters that the functions require in common.

    int n_species_0 = number of species on the first type of sites
    int n_species_1 = number of species on the second type of sites (= 0 on a 1x1 lattice)
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

    typedef double (*hmca_nn) (const int *indices, int n_species_0);

There are pre-defined functions.

    double hmca_mf_nn_1x1 (const int *indices, int n_species_0)
    double hmca_mf_nn_2x1 (const int *indices, int n_species_0)
    double hmca_pa_nn_1x1 (const int *indices, int n_species_0)
    double hmca_pa_nn_2x1 (const int *indices, int n_species_0)

### Example

Consider the Langmuir-Hinshelwood mechanism:

        O --> A
        O --> B
    A + B --> O + O

The code would look like:

    enum {O, A, B};
    int n_species_0 = 3
    int n_species_1 = 0
    int n_unimol    = 2
    int n_bimol     = 1
    int reactions[] = {O,A, O,B, A,B,O,O};
    double rates[]  = {1.0, 1.0, 10.0};

And use `hmca_mf_nn_1x1` or `hmca_pa_nn_1x1`.

## Uniform Mean-Field

    hmca_mf_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol,
        const int *reactions, const double *rates, hmca_nn nn
        )

## Half Heterogeneous Pair Approximation

    hmca_hhpa_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol, int mesh,
        const int *reactions, const double *rates, const double *weights, hmca_nn nn
        )

## Examples

To compile test_pa.c

    gcc test_pa.c hmca.c -lm
