
from matplotlib import pyplot
from numpy import linalg as la
import kwant
import numpy as np
import scipy.sparse.linalg as sla

f_nm2au = 18.89726133921252
f_eV2au = 0.03674932587122423
f_B2au=4.254382E-6
h = 1
m = 1

#############################################
############# Onside methods ################
#############################################

#     t = (h*h)/(2*m)
def onsite_e_p(sitei, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx/dx - mu
def onsite_e_m(sitei, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx/dx - mu
def onsite_h_p(sitei, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx/dx + mu
def onsite_h_m(sitei, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx/dx + mu

#############################################
############# Hopping methods ###############
#############################################

########### Hopping between +/- #############
# hopping from e+ -> e-
def hopping_e_p_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - 1j*By)
# hopping from e- -> e+
def hopping_e_m_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx + 1j*By)
# hopping from h+ -> h-
def hopping_h_p_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + 1j*By)
# hopping from h- -> h+
def hopping_h_m_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx - 1j*By)

##### Hopping +/- 1 (inside sublattice) ####
# hopping from e+ -> 1 + e+
def hopping_e(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx/dx
# hopping from h+ -> 1 + h+
def hopping_h(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx/dx

##### Hopping +/- 1 (between sublattices) ####
# hopping from e+ -> (e-) + 1
def hopping_e_p_e_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(alpha/2.0/dx)
# hopping from e+ -> (e-) - 1
def hopping_e_p_e_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (alpha/2.0/dx)
# hopping from e- -> (e+) + 1
def hopping_e_m_e_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (alpha/2.0/dx)
# hopping from e- -> (e+) - 1
def hopping_e_m_e_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(alpha/2.0/dx)
# hopping from h+ -> (h-) + 1
def hopping_h_p_h_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (alpha/2.0/dx)
# hopping from h+ -> (h-) - 1
def hopping_h_p_h_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(alpha/2.0/dx)
# hopping from h- -> (h+) + 1
def hopping_h_m_h_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(alpha/2.0/dx)
# hopping from h- -> (h+) - 1
def hopping_h_m_h_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (alpha/2.0/dx)


##### Hopping between sublattices ####
# hopping from e+ -> h+
def hopping_e_h(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta


#############################################
############### Make system #################
#############################################

def make_system(L, dx, t, mu, delta, alpha, Bx, By):
    sys = kwant.Builder()

    lat_e_p = kwant.lattice.chain(dx, name="e_p")
    lat_e_m = kwant.lattice.chain(dx, name="e_m")
    lat_h_p = kwant.lattice.chain(dx, name="h_p")
    lat_h_m = kwant.lattice.chain(dx, name="h_m")

    sys[(lat_e_p(i) for i in range(L))] = onsite_e_p
    sys[(lat_e_m(i) for i in range(L))] = onsite_e_m
    sys[(lat_h_p(i) for i in range(L))] = onsite_h_p
    sys[(lat_h_m(i) for i in range(L))] = onsite_h_m

    # hopping +/- 1 in sublattice
    sys[kwant.builder.HoppingKind((1,), lat_e_p, lat_e_p)] = hopping_e
    sys[kwant.builder.HoppingKind((-1,), lat_e_p, lat_e_p)] = hopping_e

    sys[kwant.builder.HoppingKind((1,), lat_e_m, lat_e_m)] = hopping_e
    sys[kwant.builder.HoppingKind((-1,), lat_e_m, lat_e_m)] = hopping_e

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_p)] = hopping_h
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_p)] = hopping_h

    sys[kwant.builder.HoppingKind((1,), lat_h_m, lat_h_m)] = hopping_h
    sys[kwant.builder.HoppingKind((-1,), lat_h_m, lat_h_m)] = hopping_h

    #delta
    sys[((lat_e_p(i),lat_h_m(i)) for i in range(L))] = hopping_e_h
    sys[((lat_e_m(i),lat_h_p(i)) for i in range(L))] = hopping_e_h
    sys[((lat_h_m(i),lat_e_p(i)) for i in range(L))] = hopping_e_h
    sys[((lat_h_p(i),lat_e_m(i)) for i in range(L))] = hopping_e_h

    #b
    sys[((lat_e_p(i),lat_e_m(i)) for i in range(L))] = hopping_e_p_e_m
    sys[((lat_e_m(i),lat_e_p(i)) for i in range(L))] = hopping_e_m_e_p
    sys[((lat_h_m(i),lat_h_p(i)) for i in range(L))] = hopping_h_m_h_p
    sys[((lat_h_p(i),lat_h_m(i)) for i in range(L))] = hopping_h_p_h_m

    # hopping e+ -> e- +/ 1

    sys[kwant.builder.HoppingKind((1,), lat_e_p, lat_e_m)] = hopping_e_p_e_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_e_p, lat_e_m)] = hopping_e_p_e_m_down

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_m)] = hopping_h_p_h_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_m)] = hopping_h_p_h_m_down

    sys[kwant.builder.HoppingKind((1,), lat_h_m, lat_h_p)] = hopping_h_m_h_p_up
    sys[kwant.builder.HoppingKind((-1,), lat_h_m, lat_h_p)] = hopping_h_m_h_p_down

    sys[kwant.builder.HoppingKind((1,), lat_e_m, lat_e_p)] = hopping_e_m_e_p_up
    sys[kwant.builder.HoppingKind((-1,), lat_e_m, lat_e_p)] = hopping_e_m_e_p_down

    sys = sys.finalized()
    return sys

def simulate(L, dx, t, mu, delta, alpha, Bx, By):
    sys = make_system(L, dx, t, mu, delta, alpha, Bx, By)

    ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By], sparse=False)
    e, v = la.eigh(ham_mat)
    return e

def run(L, dx, t, mu, delta, alpha, Bx, By):
    e = simulate(L, dx, t, mu, delta, alpha, Bx, By)
    return min(list(map(lambda x: abs(x), e)))

f = open('data.dat', 'w')

# L = 100
# dx = 0.1*f_nm2au
# t = 1.0/2.0/m
# mu = 10e-3*f_eV2au
# Bx = 0.0*f_B2au
# By = 0.0*f_B2au
# delta = 200e-3*f_eV2au
# alpha = 0.0*f_eV2au*f_nm2au
#
# sys = make_system(L, dx, t, mu, delta, alpha, Bx, By)
#
# ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By], sparse=False)
# e, v = la.eigh(ham_mat)
# sites = sys.sites
#
# ## Save eigenvalues
#
# eigenvalues = open('eigenvalues.dat', 'w')
# for i in range(0, len(e)):
#     eigenvalues.write("%e\t%e\n" % (i, e[i]/f_eV2au))
