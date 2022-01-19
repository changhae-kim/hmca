#/usr/bin/python3

import numpy
import hmca

def sym_id (i, j, n):
    if i < j:
        return n*i+j-i*(i+1)//2
    else:
        return n*j+i-j*(j+1)//2

def pa_closure (indices, y, n_species_0, n_species_1, params):
    n_species = n_species_0+n_species_1
    nn_species = n_species*(n_species+1)//2
    indices = numpy.ctypeslib.as_array(indices, (3,))
    y = numpy.ctypeslib.as_array(y, (nn_species,))
    ij = sym_id(indices[0], indices[1], n_species)
    jk = sym_id(indices[1], indices[2], n_species)
    ijk = nn_species*indices[1]+sym_id(indices[0], indices[2], n_species)
    yj = 0.0
    for a in range(n_species):
        yj += y[sym_id(indices[1], a, n_species)]
    return y[ij]*y[jk]/(yj+(yj == 0.0))

def pa_deriv (indices, y, n_species_0, n_species_1, params):
    n_species = n_species_0+n_species_1
    nn_species = n_species*(n_species+1)//2
    indices = numpy.ctypeslib.as_array(indices, (5,))
    y = numpy.ctypeslib.as_array(y, (nn_species,))
    ij = sym_id(indices[0], indices[1], n_species)
    jk = sym_id(indices[1], indices[2], n_species)
    ijk = nn_species*indices[1]+sym_id(indices[0], indices[2], n_species)
    lm = sym_id(indices[3], indices[4], n_species)
    yj = 0.0
    for a in range(n_species):
        yj += y[sym_id(indices[1], a, n_species)]
    return (lm == ij) * (y[jk]+(yj == 0.0 and lm == jk))/(yj+(yj == 0.0)) \
            + (lm == jk) * (y[ij]+(yj == 0.0 and lm == ij))/(yj+(yj == 0.0)) \
            - (indices[1] == indices[3] or indices[1] == indices[4]) * (y[ij]*y[jk]+(yj == 0.0 and lm == ij and lm == jk))/(yj*yj+(yj == 0.0))

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

mymf = hmca.mlmc(n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, hmca.pa_nn_2x1, pa_closure, pa_deriv)

dydt = mymf.func(y)
print('f =')
print(dydt)

dfdy = mymf.jac(y)
print('J =')
print(dfdy)

