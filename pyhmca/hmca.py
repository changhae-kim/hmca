import ctypes
import numpy
import os


# Load C Library #

libhmca = ctypes.CDLL(os.path.abspath('libhmca.so'))


# Type Def #

ctypes_hmca_nn = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_int)

ctypes_hmca_mc = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p))


# Shared Functions #

libhmca.hmca_sym_id.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int)
libhmca.hmca_sym_id.restype = ctypes.c_int

libhmca.hmca_lognorm_set.argtypes = (
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double
        )

libhmca.hmca_logexp_set.argtypes = (
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double
        )

libhmca.hmca_logsech_set.argtypes = (
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double
        )

libhmca.hmca_logpoisson2_set.argtypes = (
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double
        )

def lognorm_set (logk0, dlogk, n_unimol, n_bimol, mesh, bound):
    global libhmca
    rates = (mesh*(n_unimol+n_bimol)*ctypes.c_double)()
    weights = (mesh*ctypes.c_double)()
    libhmca.hmca_lognorm_set(
            ((n_unimol+n_bimol)*ctypes.c_double)(logk0), ((n_unimol+n_bimol)*ctypes.c_double)(dlogk),
            rates, weights,
            ctypes.c_int(n_unimol), ctypes.c_int(n_bimol), ctypes.c_int(mesh), ctypes.c_double(bound)
            )
    return rates, weights

def logexp_set (logk0, dlogk, n_unimol, n_bimol, mesh, bound):
    global libhmca
    rates = (mesh*(n_unimol+n_bimol)*ctypes.c_double)()
    weights = (mesh*ctypes.c_double)()
    libhmca.hmca_logexp_set(
            ((n_unimol+n_bimol)*ctypes.c_double)(logk0), ((n_unimol+n_bimol)*ctypes.c_double)(dlogk),
            rates, weights,
            ctypes.c_int(n_unimol), ctypes.c_int(n_bimol), ctypes.c_int(mesh), ctypes.c_double(bound)
            )
    return rates, weights

def logsech_set (logk0, dlogk, n_unimol, n_bimol, mesh, bound):
    global libhmca
    rates = (mesh*(n_unimol+n_bimol)*ctypes.c_double)()
    weights = (mesh*ctypes.c_double)()
    libhmca.hmca_logsech_set(
            ((n_unimol+n_bimol)*ctypes.c_double)(logk0), ((n_unimol+n_bimol)*ctypes.c_double)(dlogk),
            rates, weights,
            ctypes.c_int(n_unimol), ctypes.c_int(n_bimol), ctypes.c_int(mesh), ctypes.c_double(bound)
            )
    return rates, weights

def logpoisson2_set (logk0, dlogk, n_unimol, n_bimol, mesh, bound):
    global libhmca
    rates = (mesh*(n_unimol+n_bimol)*ctypes.c_double)()
    weights = (mesh*ctypes.c_double)()
    libhmca.hmca_logpoisson2_set(
            ((n_unimol+n_bimol)*ctypes.c_double)(logk0), ((n_unimol+n_bimol)*ctypes.c_double)(dlogk),
            rates, weights,
            ctypes.c_int(n_unimol), ctypes.c_int(n_bimol), ctypes.c_int(mesh), ctypes.c_double(bound)
            )
    return rates, weights


# Mean-Field Approximation #

libhmca.hmca_mf_nn_1x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_mf_nn_1x1.restype = ctypes.c_double

libhmca.hmca_mf_nn_2x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_mf_nn_2x1.restype = ctypes.c_double

libhmca.hmca_mf_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_mf_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

mf_nn_1x1 = libhmca.hmca_mf_nn_1x1

mf_nn_2x1 = libhmca.hmca_mf_nn_2x1


class mf:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        self.n_species = n_species_0+n_species_1
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.n_species*ctypes.c_double)()
        libhmca.hmca_mf_func(
                (self.n_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.n_species*self.n_species*ctypes.c_double)()
        libhmca.hmca_mf_jac(
                (self.n_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.n_species, self.n_species))


# Pair Approximation #

libhmca.hmca_pa_nn_1x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_pa_nn_1x1.restype = ctypes.c_double

libhmca.hmca_pa_nn_2x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_pa_nn_2x1.restype = ctypes.c_double

libhmca.hmca_pa_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_pa_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

pa_nn_1x1 = libhmca.hmca_pa_nn_1x1

pa_nn_2x1 = libhmca.hmca_pa_nn_2x1


class pa:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        n_species = n_species_0+n_species_1
        self.nn_species = n_species*(n_species+1)//2
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.nn_species*ctypes.c_double)()
        libhmca.hmca_pa_func(
                (self.nn_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.nn_species*self.nn_species*ctypes.c_double)()
        libhmca.hmca_pa_jac(
                (self.nn_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.nn_species, self.nn_species))


# Select Pair Approximation #

libhmca.hmca_spa_nn_1x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_spa_nn_1x1.restype = ctypes.c_double

libhmca.hmca_spa_nn_2x1.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
libhmca.hmca_spa_nn_2x1.restype = ctypes.c_double

libhmca.hmca_spa_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_spa_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

spa_nn_1x1 = libhmca.hmca_spa_nn_1x1

spa_nn_2x1 = libhmca.hmca_spa_nn_2x1


class spa:

    def __init__ (self, n_species_0, n_species_1, n_pairs, n_unimol, n_bimol, pairs, reactions, rates, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        self.n_pairs = n_pairs
        self.n_vars = n_species_0+n_species_1+n_pairs
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.pairs = pairs.copy()
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.n_vars*ctypes.c_double)()
        libhmca.hmca_spa_func(
                (self.n_vars*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1), ctypes.c_int(self.n_pairs),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                (2*self.n_pairs*ctypes.c_int)(*self.pairs),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.n_vars*self.n_vars*ctypes.c_double)()
        libhmca.hmca_spa_jac(
                (self.n_vars*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1), ctypes.c_int(self.n_pairs),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                (2*self.n_pairs*ctypes.c_int)(*self.pairs),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.n_vars, self.n_vars))


# Heterogeneous Mean-Field Approximation #

libhmca.hmca_hmf_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_hmf_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

hmf_nn_1x1 = libhmca.hmca_mf_nn_1x1

hmf_nn_2x1 = libhmca.hmca_mf_nn_2x1


class hmf:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        self.n_species = n_species_0+n_species_1
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.mesh = mesh
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.weights = weights.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.mesh*self.n_species*ctypes.c_double)()
        libhmca.hmca_hmf_func(
                (self.mesh*self.n_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.mesh*self.n_species*self.mesh*self.n_species*ctypes.c_double)()
        libhmca.hmca_hmf_jac(
                (self.mesh*self.n_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.mesh*self.n_species, self.mesh*self.n_species))


# Half Heterogeneous Pair Approximation #

libhmca.hmca_hhpa_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_hhpa_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

hhpa_nn_1x1 = libhmca.hmca_pa_nn_1x1

hhpa_nn_2x1 = libhmca.hmca_pa_nn_2x1


class hhpa:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        self.n_species = n_species_0+n_species_1
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.mesh = mesh
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.weights = weights.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.mesh*self.n_species*self.n_species*ctypes.c_double)()
        libhmca.hmca_hhpa_func(
                (self.mesh*self.n_species*self.n_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.mesh*self.n_species*self.n_species*self.mesh*self.n_species*self.n_species*ctypes.c_double)()
        libhmca.hmca_hhpa_jac(
                (self.mesh*self.n_species*self.n_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.mesh*self.n_species*self.n_species, self.mesh*self.n_species*self.n_species))


# Symmetric Heterogeneous Pair Approximation #

libhmca.hmca_shpa_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

libhmca.hmca_shpa_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn
        )

shpa_nn_1x1 = libhmca.hmca_pa_nn_1x1

shpa_nn_2x1 = libhmca.hmca_pa_nn_2x1


class shpa:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, nn):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        n_species = n_species_0+n_species_1
        self.nn_species = n_species*(n_species+1)//2
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.mesh = mesh
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.weights = weights.copy()
        self.nn = nn
        return

    def func (self, y):
        global libhmca
        dydt = (self.mesh*self.nn_species*ctypes.c_double)()
        libhmca.hmca_shpa_func(
                (self.mesh*self.nn_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.mesh*self.nn_species*self.mesh*self.nn_species*ctypes.c_double)()
        libhmca.hmca_shpa_jac(
                (self.mesh*self.nn_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol), ctypes.c_int(self.mesh),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                (self.mesh*(self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                (self.mesh*ctypes.c_double)(*self.weights),
                ctypes_hmca_nn(self.nn)
                )
        return numpy.array(dfdy).reshape((self.mesh*self.nn_species, self.mesh*self.nn_species))


# Machine Learning Moment Closure #

libhmca.hmca_mlmc_func.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn,
        ctypes_hmca_mc, ctypes_hmca_mc, ctypes.POINTER(ctypes.c_void_p)
        )

libhmca.hmca_mlmc_jac.argtypes = (
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes_hmca_nn,
        ctypes_hmca_mc, ctypes_hmca_mc, ctypes.c_void_p
        )

mlmc_nn_1x1 = libhmca.hmca_pa_nn_1x1

mlmc_nn_2x1 = libhmca.hmca_pa_nn_2x1


class mlmc:

    def __init__ (self, n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, nn, closure, deriv):
        self.n_species_0 = n_species_0
        self.n_species_1 = n_species_1
        n_species = n_species_0+n_species_1
        self.nn_species = n_species*(n_species+1)//2
        self.n_unimol = n_unimol
        self.n_bimol = n_bimol
        self.reactions = reactions.copy()
        self.rates = rates.copy()
        self.nn = nn
        self.closure = closure
        self.deriv = deriv
        return

    def func (self, y):
        global libhmca
        dydt = (self.nn_species*ctypes.c_double)()
        libhmca.hmca_mlmc_func(
                (self.nn_species*ctypes.c_double)(*y),
                dydt,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn),
                ctypes_hmca_mc(self.closure), ctypes_hmca_mc(self.deriv), ctypes.POINTER(ctypes.c_void_p)()
                )
        return numpy.array(dydt)

    def jac (self, y):
        global libhmca
        dfdy = (self.nn_species*self.nn_species*ctypes.c_double)()
        libhmca.hmca_mlmc_jac(
                (self.nn_species*ctypes.c_double)(*y),
                dfdy,
                ctypes.c_int(self.n_species_0), ctypes.c_int(self.n_species_1),
                ctypes.c_int(self.n_unimol), ctypes.c_int(self.n_bimol),
                ((2*self.n_unimol+4*self.n_bimol)*ctypes.c_int)(*self.reactions),
                ((self.n_unimol+self.n_bimol)*ctypes.c_double)(*self.rates),
                ctypes_hmca_nn(self.nn),
                ctypes_hmca_mc(self.closure), ctypes_hmca_mc(self.deriv), ctypes.POINTER(ctypes.c_void_p)()
                )
        return numpy.array(dfdy).reshape((self.nn_species, self.nn_species))


def sample_closure (indices, y, n_species_0, n_species_1, params):
    global libhmca
    n_species = n_species_0+n_species_1
    nn_species = n_species*(n_species+1)//2
    indices = numpy.ctypeslib.as_array(indices, (3,))
    y = numpy.ctypeslib.as_array(y, (nn_species,))
    ij = libhmca.hmca_sym_id(indices[0], indices[1], n_species)
    jk = libhmca.hmca_sym_id(indices[1], indices[2], n_species)
    ijk = nn_species*indices[1]+libhmca.hmca_sym_id(indices[0], indices[2], n_species)
    yj = 0.0
    for a in range(n_species):
        yj += y[libhmca.hmca_sym_id(indices[1], a, n_species)]
    return y[ij]*y[jk]/(yj+(yj == 0.0))

def sample_deriv (indices, y, n_species_0, n_species_1, params):
    global libhmca
    n_species = n_species_0+n_species_1
    nn_species = n_species*(n_species+1)//2
    indices = numpy.ctypeslib.as_array(indices, (5,))
    y = numpy.ctypeslib.as_array(y, (nn_species,))
    ij = libhmca.hmca_sym_id(indices[0], indices[1], n_species)
    jk = libhmca.hmca_sym_id(indices[1], indices[2], n_species)
    ijk = nn_species*indices[1]+libhmca.hmca_sym_id(indices[0], indices[2], n_species)
    lm = libhmca.hmca_sym_id(indices[3], indices[4], n_species)
    yj = 0.0
    for a in range(n_species):
        yj += y[libhmca.hmca_sym_id(indices[1], a, n_species)]
    return (lm == ij) * (y[jk]+(yj == 0.0 and lm == jk))/(yj+(yj == 0.0)) \
            + (lm == jk) * (y[ij]+(yj == 0.0 and lm == ij))/(yj+(yj == 0.0)) \
            - (indices[1] == indices[3] or indices[1] == indices[4]) * (y[ij]*y[jk]+(yj == 0.0 and lm == ij and lm == jk))/(yj*yj+(yj == 0.0))


if __name__ == '__main__':

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
    #n_pairs = 6
    #pairs = [0,0, 0,1, 0,2, 0,3, 1,2, 2,2]
    #mesh = 3

    rates = [float(x) for x in '''
            0.7 0.9 0.8 0.2 0.1 0.6 0.6 0.8 0.2 0.8 0.8 0.5 0.5 0.5 0.6 0.9 0.7 0.2 0.7 0.7 0.8 0.3
            '''.split()]
    #weights = [float(x) for x in '+1.00'.split()]
    y = [float(x) for x in '''
            0.00 0.33 0.00 0.61 0.30 0.00 0.93 0.19 0.90 0.00 0.21 0.21 0.00 0.09 0.02 0.30 0.44 0.00 0.97 0.25 0.34
            '''.split()]

    mymf = pa(n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, pa_nn_2x1)
    #mymf = spa(n_species_0, n_species_1, n_pairs, n_unimol, n_bimol, pairs, reactions, rates, spa_nn_2x1)
    #mymf = shpa(n_species_0, n_species_1, n_unimol, n_bimol, mesh, reactions, rates, weights, shpa_nn_2x1)
    #mymf = mlmc(n_species_0, n_species_1, n_unimol, n_bimol, reactions, rates, pa_nn_2x1, sample_closure, sample_deriv)

    dydt = mymf.func(y)
    print(dydt)

    dfdy = mymf.jac(y)
    print(dfdy)

