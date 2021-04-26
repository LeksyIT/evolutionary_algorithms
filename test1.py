import math
import matplotlib.pyplot as plt
import numpy as np

ps = 0.13
a = 1.5E5
alpha = 1.5E6
beta = 1.2E-8
k = 6E5

t_k = math.pi * 3
n = 1000
h = t_k / n
Emax = 2E6


def f1(En, P0, sig):
    pi = ps * (math.cosh((En + alpha * P0) / a) - a / (En + alpha * P0))
    return (pi - P0) / (sig * k - alpha * (pi - P0))


def RungeKutt(h):
    result = []
    P00 = 0
    E0 = 0
    sig = 0
    result.append([E0, P00])
    for i in range(int(n) - 1):
        t = h * (i + 1)
        E = Emax * math.sin(t)
        h_n = E - E0

        k1 = f1(E, P00, sig)
        k2 = f1(E + h_n / 2, P00 + k1 * h_n / 2, sig)
        k3 = f1(E + h_n / 2, P00 + k2 * h_n / 2, sig)
        k4 = f1(E + h_n, P00 + h_n * k3, sig)

        dx = (k1 + 2 * k2 + 2 * k3 + k4) / 6

        P0 = P00 + dx

        result.append([E, P0])
        P00 = P0
        sig = np.sign(h_n)
        E0 = E

    return result


E4 = RungeKutt(h)

for i in range(n ):
    E4[i][1] += beta * E4[i][0]

result = np.array(E4).T.tolist()

fig, ax = plt.subplots()
ax.plot(result[0], result[1])
plt.show()


def my_plotter(ax, data1, data2, param_dict):
    out = ax.plot(data1, data2, **param_dict)
    return out
