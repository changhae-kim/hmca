# Heterogeneous Moment Closure Approximation (HMCA)

User Guide

    hmca_mf_func (
        const double *y,
        double *dydt,
        int n_species_0, int n_species_1, int n_unimol, int n_bimol,
        const int *reactions, const double *rates, hmca_nn nn
        );


To compile test_pa.c

    gcc test_pa.c hmca.c -lm
