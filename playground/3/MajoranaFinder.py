from main import simulate
from main import f_B2au, f_nm2au, f_eV2au

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


# Sliders parameters:
mu_min = -5e-1
Bx_min = 0.0
By_min = 0.0
alpha_min = 0.0
delta_min = 0.0

mu_max = 5e-1
Bx_max = 1000e-1
By_max = 1000e-2
alpha_max = 10e-1
delta_max = 10e-1

#initial values
m = 1
h = 1

L = 100
dx = 0.1*f_nm2au
t = 1.0/2.0/m
mu = 0.1*f_eV2au
Bx = 0.2*f_B2au
By = 0.0*f_B2au
delta = 0.1*f_eV2au
alpha = 0.0*f_eV2au*f_nm2au


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.4)
X = np.linspace(0, 400, 400, endpoint=True)

plot, = plt.plot(X, simulate(L=L, dx=dx, t=t, mu=mu, Bx=Bx, By=By, delta=delta, alpha=alpha), linestyle='None', marker='o')

axcolor = 'lightgoldenrodyellow'
plt.plot([-1000, 1000], [0, 0], color='k', linestyle='-', linewidth=2)

muSlider_ax = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
BxSlider_ax = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
BySlider_ax = plt.axes([0.25, 0.20, 0.65, 0.03], axisbg=axcolor)
alphaSlider_ax = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)
deltaSlider_ax = plt.axes([0.25, 0.30, 0.65, 0.03], axisbg=axcolor)
empty_ax = plt.axes([0.25, 0.35, 0.65, 0.03], axisbg=axcolor)


muSlider = Slider(muSlider_ax, 'mu', mu_min, mu_max, valinit=(mu_min + (mu_max-mu_min)/2.0))
BxSlider = Slider(BxSlider_ax, 'Bx', Bx_min, Bx_max, valinit=(Bx_min + (Bx_max-Bx_min)/2.0))
BySlider = Slider(BySlider_ax, 'By', By_min, By_max, valinit=(By_min + (By_max-By_min)/2.0))
deltaSlider = Slider(deltaSlider_ax, 'delta', delta_min, delta_max, valinit=(delta_min + (delta_max-delta_min)/2.0))
alphaSlider = Slider(alphaSlider_ax, 'alpha', alpha_min, alpha_max, valinit=(alpha_min + (alpha_max-alpha_min)/2.0))
emptySlider = Slider(empty_ax, 'empty', 0.0, 10.0, valinit=5)

def update(val):
    mu_converted = muSlider.val*f_eV2au
    Bx_converted = BxSlider.val*f_B2au
    By_converted = BySlider.val*f_B2au
    delta_converted = deltaSlider.val*f_eV2au
    alpha_converted = alphaSlider.val*f_eV2au*f_nm2au

    data = simulate(L=L, dx=dx, t=t, mu=mu_converted, Bx=Bx_converted, By=By_converted, delta=delta_converted, alpha=alpha_converted)
    if 0.0 in data:
        print("Found!")
    plot.set_ydata(data)
    fig.canvas.draw_idle()
    plt.axis([0.0, len(data), data[0], data[-1]])
    # pyplot.draw()
deltaSlider.on_changed(update)
alphaSlider.on_changed(update)
muSlider.on_changed(update)
BxSlider.on_changed(update)
BySlider.on_changed(update)
plt.grid()
plt.show()
