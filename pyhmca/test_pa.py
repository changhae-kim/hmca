#/usr/bin/python3

import hmca

n_species_0 = 3
n_species_1 = 3
n_unimol = 4
n_bimol  = 18
reactions = [
        0,1, 3,4, 1,0, 4,3,
        0,0,2,2, 3,3,5,5, 3,0,5,2, 2,2,0,0, 5,5,3,3, 5,2,3,0,
        1,0,0,1, 4,3,3,4, 1,3,0,4, 4,0,3,1,
        2,0,0,2, 5,3,3,5, 2,3,0,5, 5,0,3,2,
        1,2,0,0, 4,5,3,3, 1,5,0,3, 4,2,3,0,
        ]

rates = [0.7, 0.9, 0.8, 0.2, 0.1, 0.6, 0.6, 0.8, 0.2, 0.8, 0.8, 0.5, 0.5, 0.5, 0.6, 0.9, 0.7, 0.2, 0.7, 0.7, 0.8, 0.3]
y = [0.00, 0.33, 0.00, 0.61, 0.30, 0.00, 0.93, 0.19, 0.90, 0.00, 0.21, 0.21, 0.00, 0.09, 0.02, 0.30, 0.44, 0.00, 0.97, 0.25, 0.34]

print('k =', rates)
print('y =', y)

mymf = hmca.pa(n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca.pa_nn_2x1)

dydt = mymf.func(y)
print('f =')
print(dydt)

dfdy = mymf.jac(y)
print('J =')
print(dfdy)

