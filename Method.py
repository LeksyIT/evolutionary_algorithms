import math
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np


def delta_func(e):
    return np.sign(2e6 * math.cos(e * h))


def cth(arg):
    return 1 / np.tanh(arg)


def p_inf(p_0, e):
    return p_s * (cth((e + alpha * p_0) / a) - a / (e + alpha * p_0))


def func_e(t):
    return E_MAX * math.sin(t)


# Обратимая часть
def func_p_e(e):
    return beta * e


# Необратимая часть
def func_p_o(p_0, t):
    return (p_inf(p_0, func_e(t)) - p_0) / (delta_func(func_e(t)) * k - alpha * (p_inf(p_0, func_e(t)) - p_0))


def rk4(p_e, t, h):
    """ Runge-Kutta 4 method """
    k1 = h * func_p_o(p_e, t)
    k2 = h * func_p_o(p_e + 0.5 * k1, t + 0.5 * h)
    k3 = h * func_p_o(p_e + 0.5 * k2, t + 0.5 * h)
    k4 = h * func_p_o(p_e + k3, t + h)
    return (k1 + 2 * k2 + 2 * k3 + k4) / 6


x = 0
p_s = 0.1  # [0.1,0.45]
beta = 1.2e-8  # [1.15e-8,1.32e-8]
a = 1e5  # [1e5,1e6]
alpha = 1e6  # [1e6,1e7]
k = 5e5  # [5e5,5e6]
E_MAX = 2e6

# t  in [0,3*math.pi]
t_start = 0
t_end = 3*math.pi
n = 100
h = (t_end - t_start) / n

p_e_points = [0]
p_0_points = [0]
p_points = [1]
t_array = [0]
e_array = [0]
test_array = [1]

for i in range(1, n):
    p_0_points.append(rk4(x, i, h))
    p_e_points.append(func_p_e(func_e(i*h)))
    t_array.append(h*i)
    e_array.append(func_e(i*h))
    p_points.append(rk4(x, i, h) + func_p_e(func_e(i*h)))
for j in range(1,n):
    test_array.append(np.sign(2e6 * math.cos(j*h)))
print(p_e_points)
fig, ax = plt.subplots()
ax.plot(e_array, p_e_points)
plt.xlabel("OE")
plt.ylabel("OP")
plt.title("Hysteresis loop")
plt.show()
