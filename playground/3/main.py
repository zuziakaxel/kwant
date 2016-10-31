
from matplotlib import pyplot
from numpy import linalg as la
import kwant
import numpy as np
import scipy.sparse.linalg as sla

f_nm2au = 18.89726133921252
f_eV2au = 0.03674932587122423
f_B2au=4.254382E-6
h = 1
m = 2

#############################################
############# Onside methods ################
#############################################

#     t = (h*h)/(2*m)
def onsite_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx - mu
def onsite_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx - mu
def onsite_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx + mu
def onsite_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx + mu

#############################################
############# Hopping methods ###############
#############################################

########### Hopping between +/- #############
# hopping from e+ -> e-
def hopping_electron_p_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) + (alpha*2*t/dx)
# hopping from e- -> e+
def hopping_electron_m_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) - (alpha*2*t/dx)
# hopping from h+ -> h-
def hopping_hole_p_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + By) + (alpha*2*t/dx)
# hopping from h- -> h+
def hopping_hole_m_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + -By) - (alpha*2*t/dx)

##### Hopping +/- 1 (inside sublattice) ####
# hopping from e+ -> 1 + e+
def hopping_electron_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx - mu
# hopping from e- -> 1 + e-
def hopping_electron_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx - mu
# hopping from h+ -> 1 + h+
def hopping_hole_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx + mu
# hopping from h- -> 1 + h-
def hopping_hole_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx + mu

#TODO:
##### Hopping +/- 1 (between sublattices) ####
# hopping from e+ -> (e-) + 1
def hopping_e_p_e_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) + (alpha*t/dx)
# hopping from e+ -> (e-) - 1
def hopping_e_p_e_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from e- -> (e+) + 1
def hopping_e_m_e_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from e- -> (e+) - 1
def hopping_e_m_e_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from h+ -> (h-) + 1
def hopping_h_p_h_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from h+ -> (h-) - 1
def hopping_h_p_h_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from h- -> (h+) + 1
def hopping_h_m_h_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
# hopping from h- -> (h+) - 1
def hopping_h_m_h_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):


##### Hopping between sublattices ####
# hopping from e+ -> h+
def hopping_electron_p_hole_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from e+ -> h-
def hopping_electron_p_hole_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from e- -> h+
def hopping_electron_m_hole_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from e- -> h-
def hopping_electron_m_hole_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from h+ -> e+
def hopping_hole_p_electron_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from h+ -> e-
def hopping_hole_p_electron_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from h- -> e+
def hopping_hole_m_electron_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from h- -> e=
def hopping_hole_m_electron_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0


#############################################
############### Make system #################
#############################################

def make_system(L, dx, Bx, By, mu, alpha):
    sys = kwant.Builder()
    # e_p -> electron up/+
    # e_m -> electron down/-

    lat_e_p = kwant.lattice.chain(dx, name="e_p")
    lat_e_m = kwant.lattice.chain(dx, name="e_m")
    lat_h_p = kwant.lattice.chain(dx, name="h_p")
    lat_h_m = kwant.lattice.chain(dx, name="h_m")

    sys[(lat_e(i) for i in range(L))] = onsite_e
    sys[(lat_h(i) for i in range(L))] = onsite_h


def make_system(L, dx, t, mu, delta):
    sys = kwant.Builder()
    lat_e = kwant.lattice.chain(dx, name="e")
    lat_h = kwant.lattice.chain(dx, name="h")

    sys[(lat_e(i) for i in range(L))] = onsite_e
    sys[(lat_h(i) for i in range(L))] = onsite_h
    # sys[lat_e.neighbors()] = hopping
    # sys[lat_h.neighbors()] = hopping
    sys[kwant.builder.HoppingKind((1,), lat_e, lat_h)] = hopping_electron_hole_plus
    # sys[kwant.builder.HoppingKind((1,), lat_h, lat_e)] = hopping_hole_electron
    sys[kwant.builder.HoppingKind((1,), lat_e, lat_e)] = hopping_electron
    sys[kwant.builder.HoppingKind((1,), lat_h, lat_h)] = hopping_hole
    # sys[kwant.builder.HoppingKind((-1,), lat_e, lat_h)] = hopping_electron_hole
    sys[kwant.builder.HoppingKind((-1,), lat_e, lat_h)] = hopping_electron_hole_minus
    sys[kwant.builder.HoppingKind((-1,), lat_e, lat_e)] = hopping_electron
    sys[kwant.builder.HoppingKind((-1,), lat_h, lat_h)] = hopping_hole

    sys = sys.finalized()
    return sys


###### CALCULATE AND PLOT ENERGY SPECTRUM
f = open('data.dat', 'w')
t = 1
L = 25
dx = 0.1*f_nm2au
delta = t
steps = 100
mu = 1
sys = make_system(L, dx, t, mu, delta)
# for i in range(steps):
#     mu = t*3*i/steps # mu/t : [0; 3]
#
#     # print("delta = %e\n mu = %e\n t = %e\nmu/t = %e" % (delta, mu, temp_t, mu/temp_t))
#
#     ham_mat = sys.hamiltonian_submatrix(args=[dx,t, mu, delta], sparse=False)
#     e, v = la.eigh(ham_mat)
#     sites= sys.sites # funckaj ktora zwraca poszczegolne wezly
#     lat = kwant.lattice.chain(dx)
#     psi= dict(zip(sites, v[:,1])) #gny chcemy stan zero
#     print(psi)
    # for i in range(L):

    #   f.write("%e %e\n" % (i*dx, np.abs(psi[lat(i)] )))
    #   f.write("%e %e\n" % (mu/t, e[i]/t))

mu = 2.5
ham_mat = sys.hamiltonian_submatrix(args=[dx,t, mu, delta], sparse=False)
e, v = la.eigh(ham_mat)
sites= sys.sites
lat_e = kwant.lattice.chain(dx, name="e")
lat_h = kwant.lattice.chain(dx, name="h")
lat = lat_e
# knocks = []
# for i in range(L):
#     knocks.append(lat_e(i))
# for i in range(L):
#     knocks.append(lat_h(i))

##### SAVE wave function of the state corresponding to the pair of Majorana modes
majorana_data = open('majorana_data.dat', 'w')
psi= dict(zip(sites, v[:,L]))
for i in range(L):
    majorana_data.write("%e %e\n" % (i,np.abs(psi[lat_e(i)].real)**2 + np.abs(psi[lat_h(i)].real)**2))

##### SAVE wave function of the first excited state
majorana_data = open('first_state.dat', 'w')
psi= dict(zip(sites, v[:,1]))
for i in range(L):
    majorana_data.write("%e %e\n" % (i,np.abs(psi[lat(i)].real)))
