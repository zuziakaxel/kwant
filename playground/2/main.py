
from matplotlib import pyplot
from numpy import linalg as la
import kwant
import numpy as np
import scipy.sparse.linalg as sla

f_nm2au = 18.89726133921252
f_eV2au = 0.03674932587122423
f_B2au=4.254382E-6

def onsite_e(sitei, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return mu
def onsite_h(sitei, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return -mu

def hopping_electron(sitei, sitej, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return -t

def hopping_hole(sitei, sitej, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return t

def hopping_hole_electron(sitei, sitej, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return delta

def hopping_electron_hole(sitei, sitej, dx, m, mu, delta):
    t = 1/2.0/m/dx/dx
    return -delta

def make_system(L, dx, m, mu, delta):
    sys = kwant.Builder()
    lat_e = kwant.lattice.chain(dx)
    lat_h = kwant.lattice.chain(dx)

    sys[(lat_e(i) for i in range(L))] = onsite_e
    sys[(lat_h(i) for i in range(L))] = onsite_h
    # sys[lat_e.neighbors()] = hopping
    # sys[lat_h.neighbors()] = hopping
    sys[kwant.builder.HoppingKind((1,), lat_e, lat_h)] = hopping_electron_hole
    # sys[kwant.builder.HoppingKind((1,), lat_h, lat_e)] = hopping_hole_electron
    sys[kwant.builder.HoppingKind((1,), lat_e, lat_e)] = hopping_electron
    sys[kwant.builder.HoppingKind((1,), lat_h, lat_h)] = hopping_hole
    # sys[kwant.builder.HoppingKind((-1,), lat_e, lat_h)] = hopping_electron_hole
    sys[kwant.builder.HoppingKind((-1,), lat_e, lat_h)] = hopping_hole_electron
    sys[kwant.builder.HoppingKind((-1,), lat_e, lat_e)] = hopping_electron
    sys[kwant.builder.HoppingKind((-1,), lat_h, lat_h)] = hopping_hole

    sys = sys.finalized()
    return sys

# def plot_spectrum(sys, muArray, delta):
#     energies = []
#     for


# L = 25
# dx = 0.1
# m = 0.067
# # mu = 0.0
# # delta = 1/2.0/m/dx/dx
#
# temp_t = 1/2.0/m/dx/dx
# mu_t = 2.0
# mu = temp_t*mu_t
# delta = temp_t
# print("delta = %e\n mu = %e\n t = %e\nmu/t = %e" % (delta, mu, temp_t, mu/temp_t))
# sys = make_system(L, dx, m, mu, delta)
# kwant.plot(sys)
# #
# ham_mat = sys.hamiltonian_submatrix(args=[dx,m, mu, delta], sparse=False)
# e, v = la.eigh(ham_mat)
# #
# for j in range(5):
#     print(j, e[j])
# #
# #
# # #poneiwaz w macierzy ham_mat polozenia wezlow nie sa po koleji - nalezy stworzyc slownik
# sites= sys.sites # funckaj ktora zwraca poszczegolne wezly
# lat = kwant.lattice.chain(dx)
# #
# psi= dict(zip(sites, v[:,1])) #gny chcemy stan zero
# #
# f = open('data.dat', 'a')
# for i in range(L):
#   # f.write("%e %e\n" % (i*dx, np.abs(psi[lat(i)].real)))
#   f.write("%e %e\n" % (mu/temp_t, e[i]/(temp_t)))



###### CALCULATE AND PLOT ENERGY SPECTRUM
f = open('data.dat', 'w')

steps = 100
for i in range(steps):
    mu_t = 3*i/steps # mu/t : [0; 3]
    L = 25
    dx = 0.1*f_nm2au
    m = 0.067

    temp_t = 1/2.0/m/dx/dx
    mu = temp_t*mu_t
    delta = temp_t
    print("delta = %e\n mu = %e\n t = %e\nmu/t = %e" % (delta, mu, temp_t, mu/temp_t))
    sys = make_system(L, dx, m, mu, delta)
    ham_mat = sys.hamiltonian_submatrix(args=[dx,m, mu, delta], sparse=False)
    e, v = la.eigh(ham_mat)
    sites= sys.sites # funckaj ktora zwraca poszczegolne wezly
    lat = kwant.lattice.chain(dx)
    psi= dict(zip(sites, v[:,1])) #gny chcemy stan zero
    for i in range(L):
      f.write("%e %e\n" % (i*dx, np.abs(psi[lat(i)] )))
    #   f.write("%e %e\n" % (mu/temp_t, e[i]/(temp_t)))
