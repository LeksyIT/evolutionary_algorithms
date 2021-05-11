import math
import math
import matplotlib.pyplot as plt
import numpy as np
import random


def mutation(i):
    if i == "0":
        return "1"
    else:
        return "0"


def pars():
    with open('test.txt') as f:
        for i in range(600):
            arr_test_k.append(float(f.readline()))


def ev_zd(k):
    ps = 0.45
    a = 0.5e5
    alpha = 0.75e6
    beta = 1.3e-8
    # k = 54.8

    t_k = math.pi * 3
    n = 600
    h = t_k / n
    Emax = 2E6

    def cth(arg):
        return 1 / np.tanh(arg)

    def f1(En, P0, delta):
        pi = ps * (cth((En + alpha * P0) / a) - a / (En + alpha * P0))
        return (pi - P0) / (delta * k * 1e4 - alpha * (pi - P0))

    def rungeKutt(h):
        result = []
        P00 = 0
        E0 = 0
        result.append([E0, P00])
        for i in range(n - 1):
            t = arr_test_i[i]
            h_n = arr_test_E[i] - E0
            delta = np.sign(h_n)
            k1 = f1(arr_test_E[i], P00, delta)
            k2 = f1(arr_test_E[i] + h_n / 2, P00 + k1 * h_n / 2, delta)
            k3 = f1(arr_test_E[i] + h_n / 2, P00 + k2 * h_n / 2, delta)
            k4 = f1(arr_test_E[i] + h_n, P00 + h_n * k3, delta)

            dx = h_n * (k1 + 2 * k2 + 2 * k3 + k4) / 6

            P0 = P00 + dx

            result.append([arr_test_E[i], P0])
            P00 = P0
            E0 = arr_test_E[i]

        return result

    E4 = rungeKutt(h)

    for i in range(n):
        E4[i][1] += beta * E4[i][0]
    return E4


def belongs(arr):
    return [int(i, 2) for i in arr if 45 <= int(i, 2) <= 90], [i for i in arr if
                                                               45 <= int(i, 2) <= 90]


def find_k(args):
    arr_k = args
    arr_k_bin = []
    for i in range(len(arr_k)):
        if arr_k[i] < 64:
            arr_k_bin.append("0" + bin(arr_k[i])[2:])
        else:
            arr_k_bin.append(bin(arr_k[i])[2:])
    new_arr_k_bin = arr_k_bin[:]
    for i in range(len(arr_k)):
        per_1 = random.randint(0, 4)
        per_2 = random.randint(0, 4)
        rand_bin = random.randint(0, 6)
        new_arr_k_bin.append(arr_k_bin[per_1][:rand_bin] + arr_k_bin[per_2][rand_bin:])
        new_arr_k_bin.append(arr_k_bin[per_2][:rand_bin] + arr_k_bin[per_1][rand_bin:])
    for i in range(5, len(new_arr_k_bin)):
        rand_mut = random.randint(0, 6)
        p = random.random()
        if p <= 0.31:
            arr = list(new_arr_k_bin[i])
            arr[rand_mut] = mutation(arr[rand_mut])
            new_arr_k_bin[i] = "".join(arr)
    arr_sum_per = []
    new_arr_k, new_arr_k_bin = belongs(new_arr_k_bin)
    for i in range(len(new_arr_k)):
        per = np.array(ev_zd(new_arr_k[i])).T.tolist()
        sum_per = 0
        for i in range(600):
            sum_per += abs(arr_test_k[i] - per[1][i])
        arr_sum_per.append(sum_per)
    arr_k.clear()
    arr_sum_per_sort = sorted(arr_sum_per)[:]
    for i in range(5):
        arr_k.append(int((new_arr_k_bin[arr_sum_per.index(arr_sum_per_sort[i])]), base=2))
    print(arr_k, arr_sum_per_sort[:5], arr_sum_per_sort[5:])
    return arr_k


def prs_E():
    with open('test_e.txt') as f:
        for i in range(600):
            arr_test_E.append(float(f.readline()))


def prs_i():
    with open('test_i.txt') as f:
        for i in range(600):
            arr_test_i.append(float(f.readline()))


def necessary_k(test_k):
    test_arr1 = np.array(ev_zd(test_k)).T.tolist()
    sumpp = 0
    for i in range(600):
        sumpp += abs(arr_test_k[i] - test_arr1[1][i])
    return sumpp


# zd_2
arr_test_E = []
arr_test_k = []
arr_test_i = []
pars()
prs_E()
prs_i()
arr_test_val = {}
test_arr = np.array(ev_zd(54.8)).T.tolist()
with open('test.txt') as f:
    for i in range(600):
        arr_test_i.append(float(f.readline()))
for j in range(45, 92):
    arr_test_val[j] = necessary_k(j)
print(arr_test_val)
arr_k = [55, 65, 75, 85, 95]
print("[90, 51, 61, 71, 81]")
for i in range(100):
    arr_k = find_k(arr_k)[:]
