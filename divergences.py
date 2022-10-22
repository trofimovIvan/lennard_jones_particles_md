import numpy as np
import matplotlib.pyplot as plt
from test_1 import create_display_list, moleculs_not_in_cube, find_nearest_display, vec_module

def get_molecs_cords_list(moment):
    particles_c = []
    with open('C:/Users/Home/PycharmProjects/project2/data/test_3/xyz cor/moment{}.xyz'.format(moment), 'r') as my_data:
        list_of_lines = my_data.readlines()
        for i in range(2, len(list_of_lines)):
            particles_c.append(np.array([float(list_of_lines[i].split()[0]), float(list_of_lines[i].split()[1]),
                                         float(list_of_lines[i].split()[2])]))
    return particles_c


def infinity_cond(part_list, num):
    my_particle = part_list[num]
    n = len(part_list)
    for i in range(n):
        if moleculs_not_in_cube(my_particle, part_list[i], 7.0):
            disp_list = create_display_list(part_list[i][0], part_list[i][1], part_list[i][2], 7.0)
            r_disp = find_nearest_display(disp_list, my_particle)
            part_list[i] = r_disp + my_particle
    return part_list


def u(particle_list, particle_num):
    R = 0.5
    dr = 0.1
    R_list = []
    n_list = []
    N = len(particle_list)
    my_partcile = particle_list[particle_num]
    while R < 3.4:
        n = 0
        for i in range(N):
            if 0 <= vec_module(my_partcile - particle_list[i]) - R <= dr:
                n += 1
        R_list.append(R)
        n_list.append(n / (4*np.pi*R**2 * dr))
        R += dr / 2
    return R_list, n_list

def u_final(part_list):
    r_list = []
    u_list = []
    n = len(part_list)
    R_lists = []
    u_lists = []
    for i in range(n):
        part_list = infinity_cond(part_list, i)
        R_list, n_list = u(part_list, i)
        R_lists.append(R_list)
        u_lists.append(n_list)
    q = len(R_lists[0])
    for i in range(q):
        r = 0
        u0 = 0
        print(i)
        for j in range(n):
            r += R_lists[j][i]
            u0 += u_lists[j][i]
        r /= n
        u0 /= n
        r_list.append(r)
        u_list.append(u0)
    return r_list, u_list



oj = get_molecs_cords_list(4886)
R_list, n_list = u_final(oj)


plt.plot(R_list, n_list)
plt.show()
