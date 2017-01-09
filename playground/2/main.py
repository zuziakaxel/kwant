
from matplotlib import pyplot
from numpy import linalg as la
import kwant
import numpy as np
import scipy.sparse.linalg as sla
import matplotlib.pyplot as plt

f_nm2au = 18.89726133921252
f_eV2au = 0.03674932587122423
f_B2au=4.254382E-6

def onsite_e(sitei, dx, t, mu, delta):
    return -mu
def onsite_h(sitei, dx, t, mu, delta):
    return mu

def hopping_electron(sitei, sitej, dx, t, mu, delta):
    return -t

def hopping_hole(sitei, sitej, dx, t, mu, delta):
    return t

def hopping_electron_hole_minus(sitei, sitej, dx, t, mu, delta):
    return -delta

def hopping_electron_hole_plus(sitei, sitej, dx, t, mu, delta):
    return delta

# def hopping_hole_electron_minus(sitei, sitej, dx, t, mu, delta):
#     return delta
#
# def hopping_electron_hole_minus(sitei, sitej, dx, t, mu, delta):
#     return -delta

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
steps = 200
mu = 1
sys = make_system(L, dx, t, mu, delta)

data_files = []
data_files_v = []
for i in range(2*L):
    data_files.append(open("dane/data%d.dat" % (i), 'w'))


for i in range(steps):
    data_files_v.append(open("dane/v/data%d.dat" % (i), 'w'))
    mu = t*4.0*i/steps # mu/t : [0; 3]

    print("delta = %e\n mu = %e\n mu/t = %e" % (delta, mu, mu/t))

    ham_mat = sys.hamiltonian_submatrix(args=[dx,t, mu, delta], sparse=False)
    e, v = la.eigh(ham_mat)
    sites= sys.sites # funckaj ktora zwraca poszczegolne wezly
    psi= dict(zip(sites, v[:,49])) #gny chcemy stan zero
    print(len(v))
    psi_m = dict(zip(sites, v[:,L]))
    lat_e = kwant.lattice.chain(dx, name="e")
    lat_h = kwant.lattice.chain(dx, name="h")
    x_s = []
    ground = []
    majorana = []
    for w in range(L):
        # data_files_v[i].write("%e")
        x_s.append(w)
        ground.append(np.abs(psi[lat_e(w)].real + psi[lat_e(w)].imag)**2 + np.abs(psi[lat_h(w)].real + psi[lat_h(w)].imag)**2)
        majorana.append(np.abs(psi_m[lat_e(w)].real)**2 + np.abs(psi_m[lat_h(w)].real)**2)
        data_files[w].write("%e %e\n" % (mu/t, e[w]/t))
        # f.write("%e %e\n" % (i*dx, np.abs(psi[lat(i)] )))
        # plt.plot(np.abs(psi[lat(i)]))
        # plt.show()
    for l in range(L):
        data_files[L+l].write("%e %e\n" % (mu/t, -e[l]/t))
    for j in range(L):
        data_files_v[i].write("%e %e\n" % (x_s[j], majorana[j]))


# mu = 2.5
# ham_mat = sys.hamiltonian_submatrix(args=[dx,t, mu, delta], sparse=False)
# e, v = la.eigh(ham_mat)
# sites= sys.sites
# lat_e = kwant.lattice.chain(dx, name="e")
# lat_h = kwant.lattice.chain(dx, name="h")
# lat = lat_e
# # knocks = []
# # for i in range(L):
# #     knocks.append(lat_e(i))
# # for i in range(L):
# #     knocks.append(lat_h(i))
#
# ##### SAVE wave function of the state corresponding to the pair of Majorana modes
# majorana_data = open('majorana_data.dat', 'w')
# psi= dict(zip(sites, v[:,L]))
# for i in range(L):
#     majorana_data.write("%e %e\n" % (i,np.abs(psi[lat_e(i)].real)**2 + np.abs(psi[lat_h(i)].real)**2))
#
# ##### SAVE wave function of the first excited state
# majorana_data = open('first_state.dat', 'w')
# psi= dict(zip(sites, v[:,1]))
# for i in range(L):
#     majorana_data.write("%e %e\n" % (i,np.abs(psi[lat(i)].real)))
