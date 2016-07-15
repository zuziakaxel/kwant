
from matplotlib import pyplot
from numpy import linalg as la
import kwant
import numpy as np
import scipy.sparse.linalg as sla

f_nm2au = 18.8976 # nanometr na
f_ev2au = 0.036749

def onsite(sitei, dx, m):
    t = 1/2.0/m/dx/dx
    return 2*t

def hopping(sitei, sitej, dx, m):
    t = 1/2.0/m/dx/dx
    return -t

def make_system(L, dx, m):
    sys = kwant.Builder()
    lat = kwant.lattice.chain(dx)
    sys[(lat(i) for i in range(L))] = onsite
    sys[lat.neighbors()] = hopping
    #sys[kwant.builder.HoppingKind((1,), lat, lat)] = hopping
    #sys[kwant.builder.HoppingKind((-1,), lat, lat)] = hopping
    sys = sys.finalized()
    return sys

L = 100
dx = 0.1*f_nm2au
m = 0.067

sys = make_system(L, dx, m)
# kwant.plot(sys)

ham_mat = sys.hamiltonian_submatrix(args=[dx,m], sparse=False)
e, v = la.eigh(ham_mat)

for j in range(5):
    print(j, e[j])


#poneiwaz w macierzy ham_mat polozenia wezlow nie sa po koleji - nalezy stworzyc slownik
sites= sys.sites # funckaj ktora zwraca poszczegolne wezly 
lat = kwant.lattice.chain(dx)

psi= dict(zip(sites, v[:,1])) #gny chcemy stan zero

f = open('data.dat', 'w')
for i in range(L):
  f.write("%e %e\n" % (i*dx, psi[lat(i)].real))
  
  

