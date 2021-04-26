import math
import math
import matplotlib.pyplot as plt
import numpy as np
import random


def ev_zd(k):
    ps = 0.45
    a = 0.5e5
    alpha = 0.75e6
    beta = 1.3e-8
    # k = 54.8

    t_k = math.pi * 3
    n = 1000
    h = t_k / n
    Emax = 2E6

    def cth(arg):
        return 1 / np.tanh(arg)

    def f1(En, P0, delta):
        pi = ps * (cth((En + alpha * P0) / a) - a / (En + alpha * P0))
        return (pi - P0) / (delta * k * 1e4 - alpha * (pi - P0))

    def RungeKutt(h):
        result = []
        P00 = 0
        E0 = 0
        result.append([E0, P00])
        for i in range(n - 1):
            t = h * (i + 1)
            E = Emax * math.sin(t)
            h_n = E - E0
            delta = np.sign(h_n)
            k1 = f1(E, P00, delta)
            k2 = f1(E + h_n / 2, P00 + k1 * h_n / 2, delta)
            k3 = f1(E + h_n / 2, P00 + k2 * h_n / 2, delta)
            k4 = f1(E + h_n, P00 + h_n * k3, delta)

            dx = h_n * (k1 + 2 * k2 + 2 * k3 + k4) / 6

            P0 = P00 + dx

            result.append([E, P0])
            P00 = P0
            E0 = E

        return result

    E4 = RungeKutt(h)

    for i in range(n):
        E4[i][1] += beta * E4[i][0]
    return E4


# fig, ax = plt.subplots()
#
# ax.plot(result[0], result[1])
# # ax.plot([i*h for i in range(n)], result[1])
# plt.show()

# zd_2
test_arr = np.array(ev_zd(54.8)).T.tolist()
test_arr1 = np.array(ev_zd(55)).T.tolist()
sumpp = 0
for i in range(1000):
    sumpp += abs(test_arr[1][i] - test_arr1[1][i])
print(sumpp)
arr_k = [35, 50, 40, 70, 60]
for i in range(100):
    arr_k_bin = []
    new_per = []
    for i in range(len(arr_k)):
        if arr_k[i] < 64:
            arr_k_bin.append("0" + bin(arr_k[i])[2:])
        else:
            arr_k_bin.append(bin(arr_k[i])[2:])


    def mut(i):
        if i == "0":
            return "1"
        else:
            return "0"


    new_per = arr_k_bin[:]
    rand_bin = random.randint(1, 5)
    rand_per_1 = random.choice(arr_k_bin)
    arr_k_bin.remove(rand_per_1)
    rand_per_1_1 = random.choice(arr_k_bin)
    arr_k_bin.remove(rand_per_1_1)
    new_per.append(rand_per_1[:rand_bin] + rand_per_1_1[rand_bin:])
    new_per.append(rand_per_1_1[:rand_bin] + rand_per_1[rand_bin:])
    rand_bin = random.randint(1, 5)
    rand_per_1 = random.choice(arr_k_bin)
    arr_k_bin.remove(rand_per_1)
    rand_per_1_1 = random.choice(arr_k_bin)
    arr_k_bin.remove(rand_per_1_1)
    new_per.append(rand_per_1[:rand_bin] + rand_per_1_1[rand_bin:])
    new_per.append(rand_per_1_1[:rand_bin] + rand_per_1[rand_bin:])
    for i in range(5, len(new_per)):
        rand_mut: int = random.randint(0, 6)
        p = random.random()
        if p <= 0.31:
            arr = list(new_per[i])
            arr[rand_mut] = mut(arr[rand_mut])
            new_per[i] = "".join(arr)
    arr_sum_per = []
    for i in range(len(new_per)):
        per = np.array(ev_zd(int(new_per[i], base=2))).T.tolist()
        sum_per = 0
        for i in range(1000):
            sum_per += abs(test_arr[1][i] - per[1][i])
        arr_sum_per.append(sum_per)
    arr_k.clear()
    arr_sum_per_sort = sorted(arr_sum_per)[:]
    for i in range(5):
        arr_k.append(int((new_per[arr_sum_per.index(arr_sum_per_sort[i])]), base=2))
    print(arr_k, arr_sum_per_sort[:5])