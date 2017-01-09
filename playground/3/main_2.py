
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
t = (h*h)/(2*m)
mu=0.5
g=2

#############################################
############# Onside methods ################
#############################################


def onsite_e_p(sitei, dx, t, mu, delta, alpha, Bx, By, Bz):
    return 2*t/dx/dx - mu + Bz
def onsite_e_m(sitei, dx, t, mu, delta, alpha, Bx, By, Bz):
    return 2*t/dx/dx - mu - Bz
def onsite_h_p(sitei, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -2*t/dx/dx + mu - Bz
def onsite_h_m(sitei, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -2*t/dx/dx + mu + Bz

#############################################
############# Hopping methods ###############
#############################################

########### Hopping between +/- #############
# hopping from e+ -> e-
def hopping_e_p_e_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (Bx - 1j*By)
# hopping from e- -> e+
def hopping_e_m_e_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (Bx + 1j*By)
# hopping from h+ -> h-
def hopping_h_p_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (-Bx + 1j*By)
# hopping from h- -> h+
def hopping_h_m_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (-Bx - 1j*By)

##### Hopping +/- 1 (inside sublattice) ####
# hopping from e+ -> 1 + e+
def hopping_e(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -t/dx/dx
# hopping from h+ -> 1 + h+
def hopping_h(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return t/dx/dx

##### Hopping +/- 1 (between sublattices) ####
# hopping from e+ -> (e-) + 1
def hopping_e_p_e_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -(alpha/2.0/dx)
# hopping from e+ -> (e-) - 1
def hopping_e_p_e_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (alpha/2.0/dx)
# hopping from e- -> (e+) + 1
def hopping_h_p_h_m_up(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return (alpha/2.0/dx)
# hopping from h+ -> (h-) - 1
def hopping_h_p_h_m_down(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -(alpha/2.0/dx)


##### Hopping between sublattices ####
# hopping from e+ -> h+
def hopping_e_p_h_m(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return -delta
def hopping_e_m_h_p(sitei, sitej, dx, t, mu, delta, alpha, Bx, By, Bz):
    return delta



#############################################
############### Make system #################
#############################################

def make_system(L, dx, t, mu, delta, alpha, Bx, By, Bz):
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

    # hopping +/- 1 in sublattice
    sys[kwant.builder.HoppingKind((1,), lat_e_p, lat_e_p)] = hopping_e
    sys[kwant.builder.HoppingKind((-1,), lat_e_p, lat_e_p)] = hopping_e

    sys[kwant.builder.HoppingKind((1,), lat_e_m, lat_e_m)] = hopping_e
    sys[kwant.builder.HoppingKind((-1,), lat_e_m, lat_e_m)] = hopping_e

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_p)] = hopping_h
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_p)] = hopping_h

    sys[kwant.builder.HoppingKind((1,), lat_h_m, lat_h_m)] = hopping_h
    sys[kwant.builder.HoppingKind((-1,), lat_h_m, lat_h_m)] = hopping_h

    sys[((lat_e_p(i),lat_e_m(i)) for i in range(L))] = hopping_e_p_e_m
    sys[((lat_h_p(i),lat_h_m(i)) for i in range(L))] = hopping_h_p_h_m


    #delta
    sys[((lat_e_p(i),lat_h_m(i)) for i in range(L))] = hopping_e_p_h_m
    sys[((lat_e_m(i),lat_h_p(i)) for i in range(L))] = hopping_e_m_h_p

    # hopping e+ -> e- +/ 1
    sys[kwant.builder.HoppingKind((1,), lat_e_p, lat_e_m)] = hopping_e_p_e_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_e_p, lat_e_m)] = hopping_e_p_e_m_down

    sys[kwant.builder.HoppingKind((1,), lat_h_p, lat_h_m)] = hopping_h_p_h_m_up
    sys[kwant.builder.HoppingKind((-1,), lat_h_p, lat_h_m)] = hopping_h_p_h_m_down



    sys = sys.finalized()
    return sys

def simulate(L, dx, t, mu, delta, alpha, Bx, By, Bz):
    sys = make_system(L, dx, t, mu, delta, alpha, Bx, By, Bz)

    ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By, Bz], sparse=False)
    e, v = la.eigh(ham_mat)
    return e



eigenvalues = open('eigenvalues.dat', 'w')

m=.046
alpha = 0.565*f_eV2au*f_nm2au
dlugosc=200*f_nm2au
energia=alpha**2*m
mu = 0.0
Bx = 0.0
By = 0.0
Bz = 0.0
t = 1.0/2.0/m
delta=50e-3*f_eV2au
L = 100
dx=dlugosc/(L+1)
print "dlugosc", dlugosc/f_nm2au
print "delta", delta/f_eV2au




sys = make_system(L, dx, t, mu, delta, alpha, Bx, By, Bz)

# REEGUALR
# pole=np.linspace(0,0.9,401)
# pole=pole*f_eV2au
#
# for i in range(len(pole)):
#   print pole[i]
#   ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, pole[i], By, Bz], sparse=False)
#   e, v = la.eigh(ham_mat)
#   #e=sla.eigsh(ham_mat, k=100, which='SM', return_eigenvectors=False)
#   eigenvalues.write("%e " % (pole[i]/f_eV2au))
#   for j in range(len(e)):
#     eigenvalues.write("%e " % (e[j]/f_eV2au))
#   eigenvalues.write("\n")
#

### 0-90 phi
import math as math
# for j in range(2, 11):
#     # print j
#     filename = "dane/b/%d.dat" % (j)
#     eigenvalues = open(filename, 'w')
#
#     r = (j/10.0)*f_eV2au
#     theta = 0.0
#     for i in range(180):
#         phi = (i/90.0)*math.pi/2.0
#
#         Bx = r*math.cos(theta)*math.cos(phi)
#         By = r*math.cos(theta)*math.sin(phi)
#         Bz = r*math.sin(theta)
#         # print "phi", phi
#         # print "Bx", Bx/f_eV2au
#         # print "By", By/f_eV2au
#         ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta, alpha, Bx, By, Bz], sparse=False)
#         e, v = la.eigh(ham_mat)
#         eigenvalues.write("%e " % ((phi*180.0)/math.pi))
#         for j in range(len(e)):
#           eigenvalues.write("%e " % (e[j]/f_eV2au))
#         eigenvalues.write("\n")


### delta

delta = np.linspace(0, 1000, 10)
delta = delta*1e-3*f_eV2au
print delta
for d in range(len(delta)):
    filename = "dane/delta/%d.dat" % (d)
    eigenvalues = open(filename, 'w')
    pole=np.linspace(0,0.9,101)
    pole=pole*f_eV2au

    for i in range(len(pole)):
      print "B=", pole[i]
      print "delta= ", delta[d]
      ham_mat = sys.hamiltonian_submatrix(args=[dx, t, mu, delta[d], alpha, pole[i], By, Bz], sparse=False)
      e, v = la.eigh(ham_mat)
      #e=sla.eigsh(ham_mat, k=100, which='SM', return_eigenvectors=False)
      eigenvalues.write("%e " % (pole[i]/f_eV2au))
      for j in range(len(e)):
        eigenvalues.write("%e " % (e[j]/f_eV2au))
      eigenvalues.write("\n")
