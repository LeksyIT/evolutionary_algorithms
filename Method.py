import math

import matplotlib.pyplot as plt
import numpy as np


# Обратимая часть
def rk4(p_e, t, h):
    """ Runge-Kutta 4 method """
    k1 = h * func(p_e, t)
    k2 = h * func(p_e + 0.5 * k1, t + 0.5 * h)
    k3 = h * func(p_e + 0.5 * k2, t + 0.5 * h)
    k4 = h * func(p_e + k3, t + h)
    return (k1 + 2 * k2 + 2 * k3 + k4) / 6


def func(p_e, t):
    e = E_MAX * math.sin(t)
    return beta * e


# Необратимая часть
def delta_func(e):
    return np.sign(2e6 * math.cos(e))

def cth(arg):
    return 1/np.tanh(arg)
def p_inf(e):
    return p_s * (cth((e + alpha * p_0)/a))

p_0 = 0
x = 0
p_s = 0.1  # [0.1,0.45]
beta = 1.15e-8  # [1.15e-8,1.32e-8]
a = 1e5  # [1e5,1e6]
alpha = 1e6  # [1e6,1e7]
k = 5e5  # [5e5,5e6]
E_MAX = 2e6

# t  in [0,3*math.pi]
t_start = 0
t_end = 3 * math.pi
n = 100
h = (t_end - t_start) / n

p_e_points = []
t_array = []

for i in range(n):
    p_e_points.append(rk4(x, i, h))
    t_array.append(h * i)
    print(t_array)
fig, ax = plt.subplots()
ax.plot(t_array, p_e_points)
plt.xlabel("Time")
plt.ylabel("Func")
plt.title("Hysteresis loop")
plt.show()
