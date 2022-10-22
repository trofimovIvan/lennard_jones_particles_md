import numpy as np
import matplotlib.pyplot as plt

def get_pressure_list(num_test, j):
    with open('C:/Users/Home/PycharmProjects/project2/data/isoterm3/test_{}/data test/pressure.txt'.format(num_test),
              'r') as pressure_file:
        preesure_list = []
        all_lines = pressure_file.readlines()
        n = len(all_lines)
        for i in range(j, n):
            preesure_list.append(float(all_lines[i].split()[0]))

    return preesure_list


def get_temp_list(num_test, j):
    with open('C:/Users/Home/PycharmProjects/project2/data/isoterm3/test_{}/data test/temprature.txt'.format(num_test),
              'r') as temp_file:
        temp_list = []
        all_lines = temp_file.readlines()
        n = len(all_lines)
        for i in range(j, n):
            temp_list.append(float(all_lines[i].split()[0]))

    return temp_list

def average(arr):
    av = 0
    n = len(arr)
    for i in range(n):
        av += arr[i]
    av /= n
    return av

def get_data_to_isoterm():
    p_list = []
    n_obr_list = []
    for i in range(1, 9):
        p_list.append(average(get_pressure_list(i, 2000)))
        n_obr_list.append(i/10)
    return p_list, n_obr_list


pressure_list, n_o_list = get_data_to_isoterm()
a, b = np.polyfit(n_o_list, pressure_list, deg=1)
X = np.arange(n_o_list[0], n_o_list[-1], 0.003)
Y = [a*x + b for x in X]
plt.errorbar(n_o_list, pressure_list, fmt='s', color='blue')
#plt.plot(X, Y, color='blue')
plt.xlabel('n', fontsize='17')
plt.ylabel('P', fontsize='17')
plt.grid()
plt.title('T = 2.0', fontsize='20')
plt.show()

