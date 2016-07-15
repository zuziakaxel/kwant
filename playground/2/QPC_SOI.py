#####################################################################################################
# Program ze stalym oddzialywaniem SOI zarowno w kontaktacj jak i nanodrucie w polu magn. (Bx,By,Bz)#
#####################################################################################################

import kwant
import tinyarray
from cmath import *
import math
from matplotlib import pyplot
import tinyarray
import numpy

f_nm2au = 18.89726133921252
f_eV2au = 0.03674932587122423
f_B2au=4.254382E-6
mu=0.5


###################################################################################################
#############################   SYSTEM i LEAD   ###################################################
###################################################################################################
def f_alp(x,y,x0,sigma,p):
  return math.exp((-((x-x0)**2/(sigma/2)**2)**p))

def f_QPC(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py):
  return Vg1*math.exp(-((x-x0qpc)**2/(lxqpc/2)**2)**px)*math.exp(-((y-y0qpc)**2/(lyqpc/2)**2)**py)+Vg2*math.exp(-((x-x0qpc)**2/(lxqpc/2)**2)**px)*math.exp(-((y+y0qpc)**2/(lyqpc/2)**2)**py)

def f_QPC_div_x(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py):
  return 0.0

def f_QPC_div_y(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py):
  return 0.0

#def f_QPC(x,y,x0qpc,omega,lqpc):
  #return math.exp(-(x-x0qpc)**2/2.0/lqpc**2)*0.5*m*omega**2*y**2

#def f_QPC_div_x(x,y,x0qpc,omega,lqpc):
  #return math.exp(-(x-x0qpc)**2/2.0/lqpc**2)*(0.5*m*omega**2*y**2)*(-(x-x0qpc)/lqpc**2)

#def f_QPC_div_y(x,y,x0qpc,omega,lqpc):
  #return math.exp(-(x-x0qpc)**2/2.0/lqpc**2)*m*omega**2*y

def onsite_up(sitei, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  t=0.5/m/dx/dx
  return 4*t+0.5*g*mu*Bz+f_QPC(xi,yi,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py)

def onsite_down(sitei, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  t=0.5/m/dx/dx
  return 4*t-0.5*g*mu*Bz+f_QPC(xi,yi,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py)

def hopping(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_up_x_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t-1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fy*ts

def hopping_up_x_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t+1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fy*ts

def hopping_up_y_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_up_y_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_down_x_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t+1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fy*ts

def hopping_down_x_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t-1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fy*ts

def hopping_down_y_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_down_y_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def onsite_up_down(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  return 0.5*g*mu*(Bx-1j*By)

def hopping_up_down_x_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return 0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fz*ts

def hopping_up_down_x_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return -0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fz*ts

def hopping_up_down_y_plus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return -1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fz*ts

def hopping_up_down_y_minus(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return 1j*0.5*(f_alp(xi,yi,x0,sigma,p)+f_alp(xj,yj,x0,sigma,p))*alfa*Fz*ts


################################################################################################
######################################## SYSTEM ################################################
################################################################################################
def make_system_SOI(W, L, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):

    lat_u = kwant.lattice.square(dx, name='up')
    lat_d = kwant.lattice.square(dx, name='down')

    sys = kwant.Builder()

    sys[(lat_u(x, y) for x in range(L) for y in range(-W+1,W))] = onsite_up
    sys[(lat_d(x, y) for x in range(L) for y in range(-W+1,W))] = onsite_down

    sys[kwant.builder.HoppingKind((1, 0), lat_u, lat_u)] = hopping_up_x_plus
    sys[kwant.builder.HoppingKind((-1, 0), lat_u, lat_u)] = hopping_up_x_minus
    sys[kwant.builder.HoppingKind((1, 0), lat_d, lat_d)] = hopping_down_x_plus
    sys[kwant.builder.HoppingKind((-1, 0), lat_d, lat_d)] = hopping_down_x_minus
    sys[kwant.builder.HoppingKind((1, 0), lat_u, lat_d)] = hopping_up_down_x_plus
    sys[kwant.builder.HoppingKind((-1, 0), lat_u, lat_d)] = hopping_up_down_x_minus

    sys[kwant.builder.HoppingKind((0, 1), lat_u, lat_u)] = hopping_up_y_plus
    sys[kwant.builder.HoppingKind((0, -1), lat_u, lat_u)] = hopping_up_y_minus
    sys[kwant.builder.HoppingKind((0, 1), lat_d, lat_d)] = hopping_down_y_plus
    sys[kwant.builder.HoppingKind((0, -1), lat_d, lat_d)] = hopping_down_y_minus
    sys[kwant.builder.HoppingKind((0, 1), lat_u, lat_d)] = hopping_up_down_y_plus
    sys[kwant.builder.HoppingKind((0, -1), lat_u, lat_d)] = hopping_up_down_y_minus

    sys[((lat_u(x, y), lat_d(x, y)) for x in range(L) for y in range(-W+1,W))] = onsite_up_down

    #Symmetry for the left leads up and down
    sym_left = kwant.TranslationalSymmetry((-dx, 0))

    ##lead up
    lead_left_up = kwant.Builder(sym_left)
    lead_left_up[(lat_u(0, y) for y in range(-W+1,W))] = onsite_up
    lead_left_up[lat_u.neighbors()] = hopping

    ##lead down
    lead_left_down = kwant.Builder(sym_left)
    lead_left_down[(lat_d(0, y) for y in range(-W+1,W))] = onsite_down
    lead_left_down[lat_d.neighbors()] = hopping


    sys.attach_lead(lead_left_up)
    sys.attach_lead(lead_left_down)
    sys.attach_lead(lead_left_up.reversed())
    sys.attach_lead(lead_left_down.reversed())
    return sys.finalized()



################################################################################################
######################################## LEADS #################################################
################################################################################################

def f_alp_lead(x,y,x0,sigma,p,xtrans):
  return math.exp((-((x-x0+xtrans)**2/(sigma/2)**2)**p))

def f_QPC_lead(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py,xtrans):
  return Vg1*math.exp(-((x-x0qpc)**2/(lxqpc/2)**2)**px)*math.exp(-((y-y0qpc)**2/(lyqpc/2)**2)**py)+Vg2*math.exp(-((x-x0qpc)**2/(lxqpc/2)**2)**px)*math.exp(-((y+y0qpc)**2/(lyqpc/2)**2)**py)

def f_QPC_div_x_lead(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py,xtrans):
  return 0.0

def f_QPC_div_y_lead(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py,xtrans):
  return 0.0

#def f_QPC_lead(x,y,x0qpc,omega,lqpc,xtrans):
  #return math.exp(-(x-x0qpc+xtrans)**2/2.0/lqpc**2)*0.5*m*omega**2*y**2

#def f_QPC_div_x_lead(x,y,x0qpc,omega,lqpc,xtrans):
  #return math.exp(-(x-x0qpc+xtrans)**2/2.0/lqpc**2)*(0.5*m*omega**2*y**2)*(-(x-x0qpc+xtrans)/lqpc**2)

#def f_QPC_div_y_lead(x,y,x0qpc,omega,lqpc,xtrans):
  #return math.exp(-(x-x0qpc+xtrans)**2/2.0/lqpc**2)*m*omega**2*y

def onsite_up_lead(sitei, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  t=0.5/m/dx/dx
  return 4*t+0.5*g*mu*Bz+f_QPC_lead(xi,yi,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py,xtrans)

def onsite_down_lead(sitei, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  t=0.5/m/dx/dx
  return 4*t-0.5*g*mu*Bz+f_QPC_lead(xi,yi,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py,xtrans)

def hopping_up_x_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t-1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fy*ts

def hopping_up_x_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t+1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fy*ts

def hopping_up_y_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_up_y_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  t=0.5/m/dx/dx
  return -t

def hopping_down_x_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t+1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fy*ts

def hopping_down_x_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t-1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fy*ts

def hopping_down_y_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def hopping_down_y_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  t=0.5/m/dx/dx
  ts=0.5/dx
  return -t

def onsite_up_down_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  return 0.5*g*mu*(Bx-1j*By)

def hopping_up_down_x_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return 0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fz*ts

def hopping_up_down_x_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return -0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fz*ts

def hopping_up_down_y_plus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return -1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fz*ts

def hopping_up_down_y_minus_lead(sitei, sitej, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):
  xi,yi=sitei.pos
  xj,yj=sitej.pos
  ts=0.5/dx
  return 1j*0.5*(f_alp_lead(xi,yi,x0,sigma,p,xtrans)+f_alp_lead(xj,yj,x0,sigma,p,xtrans))*alfa*Fz*ts


def make_lead(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans):

	lat_u = kwant.lattice.square(dx,'up')
	lat_d = kwant.lattice.square(dx,'down')
	sym_lead = kwant.TranslationalSymmetry((-dx, 0))
	lead = kwant.Builder(sym_lead)

	# build up one unit cell of the lead, and add the hoppings
	lead[(lat_u(0,y) for y in range(-W+1,W))]=onsite_up_lead
	lead[(lat_d(0,y) for y in range(-W+1,W))]=onsite_down_lead

	# hoppings in x-direction
	lead[kwant.builder.HoppingKind((1, 0), lat_u, lat_u)] = hopping_up_x_plus_lead
	lead[kwant.builder.HoppingKind((-1, 0), lat_u, lat_u)] = hopping_up_x_minus_lead
	lead[kwant.builder.HoppingKind((1, 0), lat_d, lat_d)] = hopping_down_x_plus_lead
	lead[kwant.builder.HoppingKind((-1, 0), lat_d, lat_d)] = hopping_down_x_minus_lead
	lead[kwant.builder.HoppingKind((1, 0), lat_u, lat_d)] = hopping_up_down_x_plus_lead
	lead[kwant.builder.HoppingKind((-1, 0), lat_u, lat_d)] = hopping_up_down_x_minus_lead
	# hoppings in y-directions
	lead[kwant.builder.HoppingKind((0, 1), lat_u, lat_u)] = hopping_up_y_plus_lead
	lead[kwant.builder.HoppingKind((0, -1), lat_u, lat_u)] = hopping_up_y_minus_lead
	lead[kwant.builder.HoppingKind((0, 1), lat_d, lat_d)] = hopping_down_y_plus_lead
	lead[kwant.builder.HoppingKind((0, -1), lat_d, lat_d)] = hopping_down_y_minus_lead
	lead[kwant.builder.HoppingKind((0, 1), lat_u, lat_d)] = hopping_up_down_y_plus_lead
	lead[kwant.builder.HoppingKind((0, -1), lat_u, lat_d)] = hopping_up_down_y_minus_lead

	lead[((lat_d(0, y), lat_u(0, y)) for y in range(-W+1,W))] = onsite_up_down_lead

	return  lead.finalized()

def make_lead_up(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
	lat_u = kwant.lattice.square(dx,'up')
	lat_d = kwant.lattice.square(dx,'down')
	sym_lead = kwant.TranslationalSymmetry((-dx, 0))
	lead = kwant.Builder(sym_lead)

	# build up one unit cell of the lead, and add the hoppings
	lead[(lat_u(0,y) for y in range(-W+1,W))]=onsite_up
	lead[lat_u.neighbors()] = hopping
	return  lead.finalized()

def make_lead_down(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2):
	lat_u = kwant.lattice.square(dx,'up')
	lat_d = kwant.lattice.square(dx,'down')
	sym_lead = kwant.TranslationalSymmetry((-dx, 0))
	lead = kwant.Builder(sym_lead)

	# build up one unit cell of the lead, and add the hoppings
	lead[(lat_d(0,y) for y in range(-W+1,W))]=onsite_down
	lead[lat_d.neighbors()] = hopping
	return  lead.finalized()

################################################################################################
######################################## FUNCTIONS #############################################
################################################################################################


def pasma_up_down(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, myk_limit, Nk):

  f=open("band_up.dat","w");
  lead_up=make_lead_up(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2)
  band=kwant.physics.Bands(lead_up,[dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2])

  momenta = numpy.linspace(-myk_limit*dx, myk_limit*dx, Nk)
  energies = [band(k) for k in momenta]

  for i in range(len(momenta)):
    f.write("%f " % (momenta[i]/dx*f_nm2au))
    for j in range(len(energies[0])):
      f.write("%f " % (energies[i][j]/f_eV2au))
    f.write("\n")
  f.close()

  g=open("band_down.dat","w");
  lead_down=make_lead_down(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2)
  band=kwant.physics.Bands(lead_down,[dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2])

  momenta = numpy.linspace(-myk_limit*dx, myk_limit*dx, Nk)
  energies = [band(k) for k in momenta]

  for i in range(len(momenta)):
    g.write("%f " % (momenta[i]/dx*f_nm2au))
    for j in range(len(energies[0])):
      g.write("%f " % (energies[i][j]/f_eV2au))
    g.write("\n")
  g.close()


def pasma(W, dx, m, alfa, Fz, Bx, By, Bz, myk_limit, Nk, Vg1, Vg2, xtrans):

  f=open("band"+str(int(xtrans/dx))+".dat","w");
  lead=make_lead(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans)
  band=kwant.physics.Bands(lead,[dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, xtrans])

  momenta = numpy.linspace(-myk_limit*dx, myk_limit*dx, Nk)
  energies = [band(k) for k in momenta]

  for i in range(len(momenta)):
    f.write("%f " % (momenta[i]/dx*f_nm2au))
    for j in range(len(energies[0])):
      f.write("%f " % (energies[i][j]/f_eV2au))
    f.write("\n")
  f.close()


def plot_psi(sys, energy, dx, W, L, args):
	wave_f=kwant.wave_function(sys,energy,args)
	density_lead_up=(abs(wave_f(0))**2).sum(axis=0)
	density_lead_down=(abs(wave_f(1))**2).sum(axis=0)
	sites= sys.sites

	psi_lead_up= dict(zip(sites, density_lead_up))
	psi_lead_down= dict(zip(sites, density_lead_down))

	lat_u = kwant.lattice.square(dx, name='up')
	lat_d = kwant.lattice.square(dx, name='down')

	f=open('psi'+str(energy/f_eV2au)+'.dat','w')

	for i in range(0,L):
		for j in range(-W+1,W):
			f_up=psi_lead_up[lat_u(i,j)]+psi_lead_down[lat_u(i,j)]
			f_down=psi_lead_up[lat_d(i,j)]+psi_lead_down[lat_d(i,j)]
			f_updown=f_up+f_down
			f.write("%e %e %e %e %e\n"%(i*dx/f_nm2au,j*dx/f_nm2au,f_up,f_down,f_updown))
		f.write("\n")

	f.close()


def plot_spin(sys, energy, dx, W, L, args):
	wave_f=kwant.wave_function(sys, energy, args)
	Sxup=[[0 for i in range(-W+1,W)] for j in range(L)]
	Syup=[[0 for i in range(-W+1,W)] for j in range(L)]
	Szup=[[0 for i in range(-W+1,W)] for j in range(L)]

	sites=sys.sites
	lat_u = kwant.lattice.square(dx, name='up')
	lat_d = kwant.lattice.square(dx, name='down')

	#####wpuszczamy spin up
	for k in range(len(wave_f(0))):
		wf=dict(zip(sites,wave_f(0)[k]))
		for i in range(0,L):
			for j in range(-W+1,W):
				Sxup[i][j]+=wf[lat_u(i,j)].conjugate()*wf[lat_d(i,j)]+wf[lat_d(i,j)].conjugate()*wf[lat_u(i,j)]
				Syup[i][j]+=-1j*wf[lat_u(i,j)].conjugate()*wf[lat_d(i,j)]+1j*wf[lat_d(i,j)].conjugate()*wf[lat_u(i,j)]
				Szup[i][j]+=wf[lat_u(i,j)].conjugate()*wf[lat_u(i,j)]-wf[lat_d(i,j)].conjugate()*wf[lat_d(i,j)]

		del wf

	f1=open('spin_up'+str(energy/f_eV2au)+'.dat',"w")
	for i in range(0,L):
		for j in range(-W+1,W):
			f1.write("%e %e %e %e %e %e %e %e\n" % (i*dx/f_nm2au,j*dx/f_nm2au, Sxup[i][j].real, Syup[i][j].real, Szup[i][j].real, Sxup[i][j].imag, Syup[i][j].imag, Szup[i][j].imag))
		f1.write("\n")
	f1.close()


	Sxdn=[[0 for i in range(-W+1,W)] for j in range(L)]
	Sydn=[[0 for i in range(-W+1,W)] for j in range(L)]
	Szdn=[[0 for i in range(-W+1,W)] for j in range(L)]


	#####wpuszczamy spin down
	for k in range(len(wave_f(1))):
		wf=dict(zip(sites,wave_f(1)[k]))
		for i in range(0,L):
			for j in range(-W+1,W):
				Sxdn[i][j]+=wf[lat_u(i,j)].conjugate()*wf[lat_d(i,j)]+wf[lat_d(i,j)].conjugate()*wf[lat_u(i,j)]
				Sydn[i][j]+=-1j*wf[lat_u(i,j)].conjugate()*wf[lat_d(i,j)]+1j*wf[lat_d(i,j)].conjugate()*wf[lat_u(i,j)]
				Szdn[i][j]+=wf[lat_u(i,j)].conjugate()*wf[lat_u(i,j)]-wf[lat_d(i,j)].conjugate()*wf[lat_d(i,j)]

		del wf

	f1=open('spin_down'+str(energy/f_eV2au)+'.dat',"w")
	for i in range(0,L):
		for j in range(-W+1,W):
			f1.write("%e %e %e %e %e %e %e %e\n" % (i*dx/f_nm2au,j*dx/f_nm2au, Sxdn[i][j].real, Sydn[i][j].real, Szdn[i][j].real, Sxdn[i][j].imag, Sydn[i][j].imag, Szdn[i][j].imag))
		f1.write("\n")
	f1.close()


	Sx=[[0 for i in range(-W+1,W)] for j in range(L)]
	Sy=[[0 for i in range(-W+1,W)] for j in range(L)]
	Sz=[[0 for i in range(-W+1,W)] for j in range(L)]
	for i in range(0,L):
	  for j in range(-W+1,W):
	    Sx[i][j]=Sxup[i][j]+Sxdn[i][j]
	    Sy[i][j]=Syup[i][j]+Sydn[i][j]
	    Sz[i][j]=Szup[i][j]+Szdn[i][j]

	f1=open('spin'+str(energy/f_eV2au)+'.dat',"w")
	for i in range(0,L):
		for j in range(-W+1,W):
			f1.write("%e %e %e %e %e %e %e %e\n" % (i*dx/f_nm2au,j*dx/f_nm2au, Sx[i][j].real, Sy[i][j].real, Sz[i][j].real, Sx[i][j].imag, Sy[i][j].imag, Sz[i][j].imag))
		f1.write("\n")
	f1.close()


################################################################################################
######################################## PROGRAM ###############################################
################################################################################################


m=0.0465
g=-51
dx=4*f_nm2au
W=19
L=801
alfa=0.572*f_nm2au*f_nm2au
Fz=0e3*(f_eV2au/1.0e7/f_nm2au)
Fy=200e3*(f_eV2au/1.0e7/f_nm2au)
Bx=0.0*f_B2au
By=0.0*f_B2au
Bz=0.0*f_B2au

x0=L/2*dx
sigma=300*dx
p=10

x0qpc=L/2*dx
y0qpc=W*dx
lxqpc=100*dx
lyqpc=W*dx
Vg1=0.012*f_eV2au
Vg2=0.012*f_eV2au
px=1
py=1

omega=0.0e-3*f_eV2au

f=open('f_alfa.dat',"w")
for i in range(L):
  for j in range(-W+1,W):
    x=i*dx
    y=j*dx
    f.write("%e %e %e\n" % (x/f_nm2au, y/f_nm2au, f_alp(x,y,x0,sigma,p)))
  f.write("\n")
f.close()

f=open('f_qpc.dat',"w")
for i in range(L):
  for j in range(-W+1,W):
    x=i*dx
    y=j*dx
    f.write("%e %e %e\n" % (x/f_nm2au, y/f_nm2au, f_QPC(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py)/f_eV2au))
  f.write("\n")
f.close()

f=open('f_qpc_div_x.dat',"w")
for i in range(L):
  for j in range(-W+1,W):
    x=i*dx
    y=j*dx
    f.write("%e %e %e\n" % (x/f_nm2au, y/f_nm2au, f_QPC_div_x(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py)/f_eV2au))
  f.write("\n")
f.close()

f=open('f_qpc_div_y.dat',"w")
for i in range(L):
  for j in range(-W+1,W):
    x=i*dx
    y=j*dx
    f.write("%e %e %e\n" % (x/f_nm2au, y/f_nm2au, f_QPC_div_y(x,y,x0qpc,lxqpc,y0qpc,lyqpc,Vg1,Vg2,px,py)/f_eV2au))
  f.write("\n")
f.close()



pasma_up_down(W, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2, 0.1/f_nm2au, 1000)
pasma(W, dx, m, alfa, Fz, Bx, By, Bz, 0.2/f_nm2au, 1000, Vg1, Vg2, (L/2)*dx)
pasma(W, dx, m, alfa, Fz, Bx, By, Bz, 0.2/f_nm2au, 1000, Vg1, Vg2, 0)

sys=make_system_SOI(W, L, dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2)
#kwant.plot(sys)

#Obliczenie ladunku
#energy=1.0e-3*f_eV2au
#plot_psi(sys, energy, dx, W, L, [dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2])
#plot_spin(sys, energy, dx, W, L, [dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2])


##Obliczenie konduktancji
f1=open("G.dat","w")
f2=open("T.dat","w")
dE=0.01e-3*f_eV2au
ne=500
for ie in range(ne):
   print(ie)
   energy = dE+ie*dE
   smatrix = kwant.smatrix(sys, energy, [dx, m, alfa, Fz, Bx, By, Bz, Vg1, Vg2])
   tupup=smatrix.transmission(2, 0)
   tupdown=smatrix.transmission(3, 0)
   tdownup=smatrix.transmission(2, 1)
   tdowndown=smatrix.transmission(3, 1)
   j_up=tupup+tdownup
   j_down=tdowndown+tupdown
   f2.write("%e %e %e %e %e %e\n"%(energy/f_eV2au, tupup, tupdown, tdownup, tdowndown, tupup + tupdown + tdownup + tdowndown))
   f1.write("%e %e %e %e\n"%(energy/f_eV2au, j_up, j_down, j_up+j_down))
f1.close()
f2.close()
