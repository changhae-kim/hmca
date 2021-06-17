# Heterogeneous Moment Closure Approximation (HMCA)

## Common Parameters

int n_species_0 = (1x1 lattice) number of species, (2x1 lattice) number of species on the first type of sites

int n_species_1 = (1x1 lattice) zero, (2x1 lattice) number of species on the second type of sites

int n_unimol = number of unimolecular reactions

int n_bimol = number of bimolecular reactions

int *reactions = 

double *rates =

nn =

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
