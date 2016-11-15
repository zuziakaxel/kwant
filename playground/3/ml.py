from sklearn import tree
from sklearn import datasets
from sklearn.svm import SVR

from main import simulate
from main import f_B2au, f_nm2au, f_eV2au

def calculate(mu, delta, alpha, Bx, By):
    L = 10
    dx = 0.1*f_nm2au
    t = 1/2.0
    mu_converted = mu*f_eV2au
    Bx_converted = Bx*f_B2au
    By_converted = By*f_B2au
    delta_converted = delta*f_eV2au
    alpha_converted = alpha*f_eV2au*f_nm2au
    data = simulate(L=L, dx=dx, t=t, mu=mu_converted, Bx=Bx_converted, By=By_converted, delta=delta_converted, alpha=alpha_converted)
    e = min(map(lambda x: abs(x), data))
    return e



from random import uniform
import numpy as np
features = []
labels = []
f = open("test-data.dat", mode='a')
for i in range(1000):
    if i%10 == 0:
        print(i)
    mu = uniform(0, 1)
    delta = uniform(0, 1)
    alpha = uniform(0, 1)
    Bx = uniform(0, 1)
    By = uniform(0, 1)
    e = calculate(mu, delta, alpha, Bx, By)
#     features.append([mu, delta, alpha, Bx, By])
#     labels.append(e)
    f.write("%f\t %f\t%f\t%f\t%f\t%f\n" %(e, mu, delta, alpha, Bx, By))
    # print("%f,%f, %f, %f, %f, %f" %(e, mu, delta, alpha, Bx, By))
lines = list(open('test-data.dat', mode='r').read().split('\n'))
lines.pop(-1)
data = list(map(lambda line: list(map(float, line.split('\t'))), lines))
features = list(map(lambda x: [x[1], x[2], x[3], x[4], x[5]], data))
labels = list(map(lambda x: x[0], data))
features = np.array(features)
clf = SVR(kernel='rbf', C=1e3, gamma=0.1)
clf.fit(features, labels)
print(calculate(0.1, 0.1, 0.1, 0.1, 0.1))
print(clf.predict(features))
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import GradientBoostingRegressor
print(MultiOutputRegressor(GradientBoostingRegressor(random_state=0)).fit(features, labels).predict(features))
