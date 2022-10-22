import numpy as np
import matplotlib.pyplot as plt
from test_1 import *

def get_particles_velocity(moment):
    particles_velocity = []
    with open('C:/Users/Home/PycharmProjects/project2/data/isoterm4/test_1'
                       '/xyz cor/moment{}.xyz'.format(moment), 'r') as moment_file:
        list_of_lines = moment_file.readlines()
        for i in range(2, len(list_of_lines)):
            particles_velocity.append(np.array([float(list_of_lines[i].split()[3]), float(list_of_lines[i].split()[4]),
                                       float(list_of_lines[i].split()[5])]))
    return particles_velocity

def get_v_list():
    main_list = []
    for moment in range(1500, 2000):
        v_list_p = get_particles_velocity(moment)
        for i in range(len(v_list_p)):
            v_list_p[i] = np.sqrt(v_list_p[i][0] ** 2 + v_list_p[i][1] ** 2 + v_list_p[i][2] ** 2)
        main_list.append(v_list_p)
    all = []
    for lst in main_list:
        all = all + lst
    return all

v_list = get_v_list()

plt.hist(v_list, bins=50, density=True)
plt.show()