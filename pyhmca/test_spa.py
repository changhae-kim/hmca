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
n_pairs = 6
pairs = [0,0, 0,1, 0,2, 0,3, 1,2, 2,2]

rates = [0.9, 1.0, 0.2, 0.7, 0.4, 0.3, 0.1, 1.0, 0.9, 0.1, 0.5, 1.0, 0.2, 0.4, 0.3, 0.7, 0.6, 0.5, 0.5, 0.5, 0.9, 0.2]
y = [0.00, 0.05, 0.14, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.97, 0.23]

print('k =', rates)
print('y =', y)

mymf = hmca.spa(n_species_0, n_species_1, n_pairs, n_unimol, n_bimol, pairs, reactions, rates, hmca.spa_nn_2x1)

dydt = mymf.func(y)
print('f =')
print(dydt)

dfdy = mymf.jac(y)
print('J =')
print(dfdy)

