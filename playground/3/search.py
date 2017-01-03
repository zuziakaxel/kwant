from main import simulate
from main import f_B2au, f_nm2au, f_eV2au
import numpy as np


m=1
alpha = 0.0
delta = 0.1
L = 100
dx = 0.1*f_nm2au
t = 1.0/2.0/m
mu = 0.0*f_eV2au
Bx = 0.0*f_B2au
By = 0.0*f_B2au

i = 0
data_x = []
data_a = []
data_e0 = []
data_e1 = []
data_e2 = []
data_e3 = []
data_e4 = []
data_e5 = []

bx_s = [(float(i)/10000000000.0) for i in range(1000000500,1000005000)]
for bx in bx_s:

    i += 1
    # dx = L_so/L
    if i%10 == 0:
        print(i)
    es = simulate(L=L, dx=dx, t=t, mu=mu, Bx=0.0, By=bx, delta=delta, alpha=alpha)
    # min(e, key=func)
    # print(e)
    # print(min(map(lambda x: abs(x), e)))
    # print("%9.f" %(min(e)))
    # data.append((min(map(lambda x: abs(x), e)), alpha, delta))
    # data_x.append(bx/E_so)
    data_x.append(bx)
    data_a.append(alpha)
    # data_y.append(by)
    # e = list(map(lambda x: abs(x), es))

    # e.sort()
    # eigen = min(map(lambda x: abs(x), es))
    # print(len(es))
    # data_e.append(eigen)
    min((abs(x), x) for x in es)[1]
    data_e0.append(min(es, key=abs))
    # data_e0.append(-es[(L*2)])
    data_e1.append(es[(L*2) + 1])
    data_e2.append(-es[(L*2) - 1])
    data_e3.append(es[(L*2) + 2])
    data_e4.append(-es[(L*2) -2])
    data_e5.append(es[(L*2) - 3])
    # print(min(map(lambda x: abs(x), e)))
    for e in es:
        if abs(e) <= 1e-11:
            print("Znaleziono: %.9f" % e)
# with open('data.dat', mode='w') as f:
#     for e,a,d in data:
#         f.write("%.9f\t%.9f\t%.9f\n" % (e,a,d))

# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt
# import numpy as np
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X, Y, Z = axes3d.get_test_data(0.05)
# ax.plot_wireframe(data_a, data_x, data_e, rstride=10, cstride=10)
#
# ax.set_xlabel('Alpha')
# ax.set_ylabel('Bx')
# ax.set_zlabel('E')
# plt.show()
import matplotlib.pyplot as plt

# plt.axis([0, 6, 0, 20])

# fig, ax = plt.subplots()

plt.scatter(data_x, data_e0, label='e0')
# plt.plot(data_x, data_e1, label='e1')
# plt.plot(data_x, data_e2, label='e2')
# plt.plot(data_x, data_e3, label='e3')
# plt.plot(data_x, data_e4, label='e4')
# plt.plot(data_x, data_e5, label='e5')
plt.grid()
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# ax.set_xlabel('Bx')
# ax.set_ylabel('E')
plt.show()
