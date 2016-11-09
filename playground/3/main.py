
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
def onsite_e_p(sitei, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx - mu
def onsite_e_m(sitei, dx, t, mu, delta, alpha, Bx, By):
    return -2*t/dx - mu
def onsite_h_p(sitei, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx + mu
def onsite_h_m(sitei, dx, t, mu, delta, alpha, Bx, By):
    return 2*t/dx + mu

#############################################
############# Hopping methods ###############
#############################################

########### Hopping between +/- #############
# hopping from e+ -> e-
def hopping_e_p_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) + (alpha*2*t/dx)
# hopping from e- -> e+
def hopping_e_m_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) - (alpha*2*t/dx)
# hopping from h+ -> h-
def hopping_h_p_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + By) + (alpha*2*t/dx)
# hopping from h- -> h+
def hopping_h_m_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + -By) - (alpha*2*t/dx)

##### Hopping +/- 1 (inside sublattice) ####
# hopping from e+ -> 1 + e+
def hopping_e_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx - mu
def hopping_e_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx - mu
# hopping from e- -> 1 + e-
def hopping_e_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx - mu
def hopping_e_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx - mu
# hopping from h+ -> 1 + h+
def hopping_h_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx + mu
def hopping_h_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx + mu
# hopping from h- -> 1 + h-
def hopping_h_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return t/dx + mu
def hopping_h_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -t/dx + mu

#TODO:
##### Hopping +/- 1 (between sublattices) ####
# hopping from e+ -> (e-) + 1
def hopping_e_p_e_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) - (alpha*t/dx)
# hopping from e+ -> (e-) - 1
def hopping_e_p_e_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx - By) + (alpha*t/dx)
# hopping from e- -> (e+) + 1
def hopping_e_m_e_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx + By) + (alpha*t/dx)
# hopping from e- -> (e+) - 1
def hopping_e_m_e_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (Bx + By) - (alpha*t/dx)
# hopping from h+ -> (h-) + 1
def hopping_h_p_h_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + By) + (alpha*t/dx)
# hopping from h+ -> (h-) - 1
def hopping_h_p_h_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return (-Bx + By) - (alpha*t/dx)
# hopping from h- -> (h+) + 1
def hopping_h_m_h_p_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(Bx + By) + (alpha*t/dx)
# hopping from h- -> (h+) - 1
def hopping_h_m_h_p_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return -(Bx + By) - (alpha*t/dx)


##### Hopping between sublattices ####
# hopping from e+ -> h+
def hopping_e_p_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from e+ -> h-
def hopping_e_p_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from e- -> h+
def hopping_e_m_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from e- -> h-
def hopping_e_m_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from h+ -> e+
def hopping_h_p_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0
# hopping from h+ -> e-
def hopping_h_p_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from h- -> e+
def hopping_h_m_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return delta
# hopping from h- -> e-
def hopping_h_m_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By):
    return 0


#############################################
############### Make system #################
#############################################

def make_system(L, dx, t, mu, delta, alpha, Bx, By):
    sys = kwant.Builder()
    # e_p -> electron up/+
    # e_m -> electron down/-

    lat_e_p = kwant.lattice.chain(dx, name="e_p")
    lat_e_m = kwant.lattice.chain(dx, name="e_m")
    lat_h_p = kwant.lattice.chain(dx, name="h_p")
    lat_h_m = kwant.lattice.chain(dx, name="h_m")

    sys[(lat_e_p(i) for i in range(L))] = onsite_e_p
    sys[(lat_e_m(i) for i in range(L))] = onsite_e_m
    sys[(lat_h_p(i) for i in range(L))] = onsite_h_p
    sys[(lat_h_m(i) for i in range(L))] = onsite_h_m

    # hopping +/- in sublattice
    sys[kwant.builder.HoppingKind((1,), lat_e_p, lat_e_p)] = hopping_e_p_up
    sys[kwant.builder.HoppingKind((-1,), lat_e_p, lat_e_p)] = hopping_e_p_down

    sys[kwant.builder.HoppingKind((1,), lat_e_m, lat_e_m)] = hopping_e_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_e_m, lat_e_m)] = hopping_e_m_down

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_p)] = hopping_h_p_up
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_p)] = hopping_h_p_down

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_p)] = hopping_h_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_p)] = hopping_h_m_down

    # hopping between sublattices
    sys[kwant.builder.HoppingKind((0,), lat_e_p, lat_e_m)] = hopping_e_p_e_m
    sys[kwant.builder.HoppingKind((0,), lat_e_m, lat_e_p)] = hopping_e_m_e_p

    sys[kwant.builder.HoppingKind((0,), lat_h_p, lat_h_m)] = hopping_h_p_h_m
    sys[kwant.builder.HoppingKind((0,), lat_h_m, lat_h_p)] = hopping_h_m_h_p

    sys[kwant.builder.HoppingKind((0,), lat_e_p, lat_h_p)] = hopping_e_p_h_p
    sys[kwant.builder.HoppingKind((0,), lat_h_p, lat_e_p)] = hopping_h_p_e_p

    sys[kwant.builder.HoppingKind((0,), lat_e_m, lat_h_m)] = hopping_e_m_h_m
    sys[kwant.builder.HoppingKind((0,), lat_h_m, lat_e_m)] = hopping_h_m_e_m

    sys[kwant.builder.HoppingKind((0,), lat_e_p, lat_h_m)] = hopping_e_p_h_m
    sys[kwant.builder.HoppingKind((0,), lat_h_m, lat_e_p)] = hopping_h_m_e_p

    sys[kwant.builder.HoppingKind((0,), lat_e_m, lat_h_p)] = hopping_e_m_h_p
    sys[kwant.builder.HoppingKind((0,), lat_h_p, lat_e_m)] = hopping_h_p_e_m

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

###### CALCULATE AND PLOT ENERGY SPECTRUM
f = open('data.dat', 'w')

L = 10
dx = 0.1
t = 1
mu = 0.1
Bx = 0.2
By = 0.0 + 0.j
delta = 0.1
alpha = 0.4

sys = make_system(L, dx, t, mu, delta, alpha, Bx, By)
# kwant.plot(sys)
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

# mu = 2.5
ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By], sparse=False)
e, v = la.eigh(ham_mat)
sites= sys.sites
lat_e_p = kwant.lattice.chain(dx, name="e_p")
lat_e_m = kwant.lattice.chain(dx, name="e_m")
lat_h_p = kwant.lattice.chain(dx, name="h_p")
lat_h_m = kwant.lattice.chain(dx, name="h_m")
# lat = lat_e
# knocks = []
# for i in range(L):
#     knocks.append(lat_e(i))
# for i in range(L):
#     knocks.append(lat_h(i))

##### SAVE wave function of the state corresponding to the pair of Majorana modes
majorana_data = open('majorana_data.dat', 'w')
psi= dict(zip(sites, v[:,L]))
for i in range(L):
    psi_e_p = np.abs(psi[lat_e_p(i)].real)**2
    psi_e_m = np.abs(psi[lat_e_m(i)].real)**2
    psi_h_p = np.abs(psi[lat_h_p(i)].real)**2
    psi_h_m = np.abs(psi[lat_h_m(i)].real)**2
    majorana_data.write("%e %e\n" % (i, psi_e_p + psi_e_m + psi_h_p + psi_h_m))
## Save eigenvalues

eigenvalues = open('eigenvalues.dat', 'w')
for i in range(0, len(e)):
# for eigenvalue in e:
    eigenvalues.write("%e\t%e\n" % (i, e[i]))
#
# values = []
# for i in range(-10,10):
#     di = i/10.0
#     t = di
#     sys = make_system(L, dx, t, mu, delta, alpha, Bx, By)
#     ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By], sparse=False)
#     e, v = la.eigh(ham_mat)
#     # sort eigenvalues with lattice points
#     e_sorted = dict(zip(sites, e[:]))
#     e_sorted = e
#     # TODO:
#     values.append([e_sorted[0], e_sorted[1], e_sorted[2], e_sorted[3] ])
# #
# # index = 0
# # eigenvalues = open('eigenvalues.dat', 'w')
# first = list(map(lambda x: x[0], values))
# second = list(map(lambda x: x[1], values))
# third = list(map(lambda x: x[2], values))
# fourth = list(map(lambda x: x[3], values))
# eigenvalues = open('eigenvalues.dat', 'w')
# for i in range(-10,10):
#     eigenvalues.write("%e\t" %(i))
#     for j in [first, second, third, fourth]:
#
#         di = i/10.0
#         eigenvalues.write("%e\t" % (j[i]))
#     eigenvalues.write("\n")
    # eigenvalues.write("%e\t%e\t%e\t%e\t%e\n" %(di, values[i][0], values[i][1], values[i][2], values[i][3]))

##### SAVE wave function of the first excited state
# majorana_data = open('first_state.dat', 'w')
# psi= dict(zip(sites, v[:,1]))
# for i in range(L):
#     majorana_data.write("%e %e\n" % (i,np.abs(psi[lat(i)].real)))
