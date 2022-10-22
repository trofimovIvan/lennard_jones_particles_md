import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def get_pressure_list(isoterm_num, k):
    pressure_list = []
    time_list = []
    with open('C:/Users/Home/PycharmProjects/project2/data/isoterm{}/'
              'test_{}/data test/pressure.txt'.format(isoterm_num, k), 'r') as data_file:
        lines = data_file.readlines()
        n = len(lines)
        time = 0
        for i in range(n):
            if i > 2000:
                pressure_list.append(float(lines[i]))
                time_list.append(time)
                time += 0.001

    return time_list, pressure_list

def average_of_list(array):
    av = 0
    n = len(array)
    for i in range(n):
        av += array[i]
    av /= n
    return av

def get_energy_list():
    energy_list = []
    time_list = []
    with open('C:/Users/Home/PycharmProjects/project2/data/isoterm4/test_2/energy.txt', 'r') as data_file:
        lines = data_file.readlines()
        n = len(lines)
        time = 0
        for i in range(n):
            if i > 2000:
                energy_list.append(float(lines[i]))
                time_list.append(time)
                time += 0.001
    return time_list, energy_list

def draw_isoterm(n, is_num):
    pressure_list = []
    n_list = []
    m = 0.45
    for i in range(1, n + 1):
        pressure_list.append(average_of_list(get_pressure_list(4, i)[1]))
        if i <= 8:
            n_list.append(i / 10)
        else:
            n_list.append(m)
            m += 0.1
    n_list.sort()
    pressure_list.sort()
    popt, pcov = (curve_fit(func, n_list, pressure_list))
    print(popt)
    y = [func(n_list[i], popt[0], popt[1], popt[2]) for i in range(len(n_list))]
    print(np.sqrt(np.diag(pcov)))
    perr = np.sqrt(np.diag(pcov))
    T_k = 8*popt[0]/(27*popt[1])
    sigma_T_k = T_k * np.sqrt((perr[0] / popt[0])**2 + (perr[1] / popt[1])**2)
    print(T_k, sigma_T_k)
    return n_list, pressure_list, y, popt[2]



def draw_isoterm_2(k):
    pressure_list = []
    n_list = []
    n = 0.05
    for i in range(1, k + 1):
        pressure_list.append(average_of_list(get_pressure_list(5, i)[1]))
        n_list.append(n)
        n += 0.05
    pressure_list.append(5.34)
    n_list.append(0.65)
    pressure_list.append(7.81)
    n_list.append(0.7)
    n_list_2, pressure_list_2, y_2, T_2 = draw_isoterm(12, 4)
    popt, pcov = (curve_fit(func, n_list, pressure_list))
    print(popt)
    y = [func(n_list[i], popt[0], popt[1], popt[2]) for i in range(len(n_list))]
    #print(np.sqrt(np.diag(pcov)))
    perr = np.sqrt(np.diag(pcov))
    T_k = 8 * popt[0] / (27 * popt[1])
    sigma_T_k = T_k * np.sqrt((perr[0] / popt[0]) ** 2 + (perr[1] / popt[1]) ** 2)
    print(T_k, sigma_T_k)

    plt.errorbar(n_list, pressure_list, fmt='s')
    plt.plot(n_list, y)
    plt.errorbar(n_list_2, pressure_list_2, fmt='s')
    plt.plot(n_list_2, y_2)
    plt.xlabel('n', fontsize='17')
    plt.ylabel('P', fontsize='17')
    plt.title('T_1 = {}, T_2 = {}'.format(popt[2], T_2))
    # plt.title(r'$T_k$ = {}'.format(round(8*popt[0]/(27*popt[1]), 3)))
    plt.grid()
    plt.show()



def func(n, a, b, T):
    return n * T / (1 - n*b) - a*n**2

def sigma_p(p_list):
    n = len(p_list)
    sigma = 0
    av_p = average_of_list(p_list)
    for i in range(n):
        sigma += (p_list[i] - av_p)**2
    return np.sqrt(sigma / (n*(n - 1)))

def result_func():
    pressure_list_1 = []
    pressure_list_2 = []
    n_list_1 = []
    n_list_2 = []
    n = 0.1
    m = 0.45
    for i in range(1, 13):
        pressure_list_1.append(average_of_list(get_pressure_list(4, i)[1]))
        if i <= 8:
            n_list_1.append(i / 10)
        else:
            n_list_1.append(m)
            m += 0.1

    n = 0.00
    for i in range(1, 15):
        n += 0.05
        pressure_list_2.append(average_of_list(get_pressure_list(5, i)[1]))
        n_list_2.append(n)


    n_list_1.sort()
    pressure_list_1.sort()
    print(pressure_list_1)
    print(n_list_1)
    print(pressure_list_2)
    print(n_list_2)
    popt, pcov = curve_fit(func, n_list_1, pressure_list_1)
    print(popt)
    print(np.sqrt(np.diag(pcov)))
    popt2, pcov2 = curve_fit(func, n_list_2, pressure_list_2)
    print(popt2)
    print(np.sqrt(np.diag(pcov2)))

    print(8*popt[0]/ (27*popt[1]), 8*popt2[0] / (27*popt2[1]))
    n_fit_list_1 = np.arange(n_list_1[0], n_list_1[-1], 0.003)
    y_fit_list_1 = [func(n_fit_list_1[i], popt[0], popt[1], popt[2]) for i in range(len(n_fit_list_1))]
    n_fit_list_2 = np.arange(n_list_2[0], n_list_2[-1], 0.003)
    y_fit_list_2 = [func(n_fit_list_2[i], popt2[0], popt2[1], popt2[2]) for i in range(len(n_fit_list_2))]
    plt.errorbar(n_list_1, pressure_list_1, fmt='s', color='blue', label='T = 2.0 data')
    plt.plot(n_fit_list_1, y_fit_list_1, color='orange', label='T = 2.0 fit')

    plt.errorbar(n_list_2, pressure_list_2, fmt='s', color='red', label='T = 3.0 data')
    plt.plot(n_fit_list_2, y_fit_list_2, color='green', label='T = 3.0 fit')
    plt.grid()
    plt.legend(fontsize='13')
    plt.xlabel('n', fontsize='17')
    plt.ylabel('P', fontsize='17')
    plt.show()


"""time_list, pressure_list = get_pressure_list(6, 1)
plt.plot(time_list, pressure_list)
plt.show()
print(average_of_list(get_pressure_list(6, 1)[1]))"""

result_func()

