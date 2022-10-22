import numpy as np
import matplotlib.pyplot as plt
import random

""" Here functions to make simulation lennard-jones particle simulation
        using langevin NVT ensemble
""""

#constants

T = 1.9
gamma = 1
s = np.sqrt(2*T*gamma)


def get_cm_v(v_list):
    n = len(v_list)
    v_cm = np.array([0, 0 ,0])
    for i in range(n):
        v_cm = v_cm + v_list[i]

    v_cm = v_cm/n
    return v_cm

def get_norm_molec_v(v_list):
    n = len(v_list)
    v_cm = get_cm_v(v_list)
    for i in range(n):
        v_list[i] = v_list[i] - v_cm
    return v_list


def vec_module(array):
    return np.sqrt(array[0] ** 2 + array[1] ** 2 + array[2] ** 2)

def spawn_moleculs(part_list, l, numbers_particles, sigma):
    for i in range(int(numbers_particles)):
        x_cor = random.uniform(0, l)
        y_cor = random.uniform(0, l)
        z_cor = random.uniform(0, l)

        particle = np.array([x_cor, y_cor, z_cor])
        j = 0
        while j < len(part_list):
            vec_between = particle - part_list[j]
            displays = create_display_list(particle[0], particle[1], particle[2], l)
            nearest_disp_dist = find_nearest_display(displays, part_list[j])
            if vec_module(vec_between) <= sigma - 0.2 or vec_module(nearest_disp_dist) <= sigma - 0.2:
                x_cor = random.uniform(0, l)
                y_cor = random.uniform(0, l)
                z_cor = random.uniform(0, l)
                particle = np.array([x_cor, y_cor, z_cor])
                j = -1
            j += 1
        part_list.append(particle)
        print(i)
    return part_list

def get_velocity_to_particles(part_list, temprature):
    velocity_list = []
    for i in range(len(part_list)):
        v_x = np.sqrt(temprature) * random.choice([-1, 1])
        v_y = np.sqrt(temprature) * random.choice([-1, 1])
        v_z = np.sqrt(temprature) * random.choice([-1, 1])

        velocity_particle = np.array([v_x, v_y, v_z])
        velocity_list.append(velocity_particle)
    return velocity_list

def create_display_list(x, y, z, l):
    display_list = []
    display_list.append(np.array([x + l, y, z]))
    display_list.append(np.array([x + l, y, z + l]))
    display_list.append(np.array([x + l, y + l, z]))
    display_list.append(np.array([x + l, y + l, z + l]))
    display_list.append(np.array([x + l, y - l, z + l]))
    display_list.append(np.array([x + l, y - l, z]))
    display_list.append(np.array([x + l, y - l, z - l]))
    display_list.append(np.array([x + l, y, z - l]))
    display_list.append(np.array([x + l, y + l, z - l]))

    display_list.append(np.array([x - l, y, z]))
    display_list.append(np.array([x - l, y, z + l]))
    display_list.append(np.array([x - l, y + l, z]))
    display_list.append(np.array([x - l, y + l, z + l]))
    display_list.append(np.array([x - l, y - l, z + l]))
    display_list.append(np.array([x - l, y - l, z]))
    display_list.append(np.array([x - l, y - l, z - l]))
    display_list.append(np.array([x - l, y, z - l]))
    display_list.append(np.array([x - l, y + l, z - l]))

    display_list.append(np.array([x, y, z]))
    display_list.append(np.array([x, y, z + l]))
    display_list.append(np.array([x, y + l, z]))
    display_list.append(np.array([x, y + l, z + l]))
    display_list.append(np.array([x, y - l, z + l]))
    display_list.append(np.array([x, y - l, z]))
    display_list.append(np.array([x, y - l, z - l]))
    display_list.append(np.array([x, y, z - l]))
    display_list.append(np.array([x, y + l, z - l]))
    return display_list

def moleculs_not_in_cube(part_1, part_2, l):
    return abs(part_1[0] - part_2[0]) > l / 2 or abs(part_1[1] - part_2[1]) > l / 2 or abs(
        part_1[2] - part_2[2]) > l / 2


def find_nearest_display(disp_list, particle):
    r_disp = disp_list[1]
    min_dist = vec_module(r_disp - particle)
    for i in range(len(disp_list)):
        r_cur = disp_list[i]
        r = vec_module(r_cur - particle)
        if r <= min_dist:
            min_dist = r
            r_disp = r_cur
    return r_disp - particle


def calc_force_energy(part_list, v_list, r_cut, dt, l):
    kinetic = 0
    potent = 0

    n = len(part_list)
    force_list = [np.array([0, 0, 0]) for i in range(n)]
    for i in range(n):
        kinetic += vec_module(v_list[i])**2 / 2
    temprature = kinetic*2 / (3*n)
    vir = 0
    for i in range(n):
        for j in range(i+1,n):
            if not moleculs_not_in_cube(part_list[i], part_list[j], l):
                r = vec_module(part_list[i] - part_list[j])
                if r > r_cut:
                    continue
                potent += 4 * (r ** -12 - r ** -6)
                dudr = -24 * (-2 * r ** -13 + r ** -7)
                force_list[i] = force_list[i] + dudr * (part_list[i] - part_list[j]) / r
                force_list[j] = force_list[j] - dudr * (part_list[i] - part_list[j]) / r
                vir += np.dot(dudr * (part_list[i] - part_list[j]) / r, part_list[i] - part_list[j])
            else:
                disp_list = create_display_list(part_list[j][0], part_list[j][1], part_list[j][2], l)
                r_vec = find_nearest_display(disp_list, part_list[i])
                r = vec_module(r_vec)
                if r > r_cut:
                    continue
                potent += 4 * (r ** -12 - r ** -6)
                dudr = -24 * (-2 * r ** -13 + r ** -7)
                force_list[i] = force_list[i] - dudr * (r_vec) / r
                force_list[j] = force_list[j] + dudr * (r_vec) / r
                vir += np.dot(dudr * r_vec / r, r_vec)

    return  kinetic, potent, force_list, vir

def R(force_list, dt):
    n = len(force_list)
    for i in range(n):
        force_list[i][0] += random_force(T, dt)
        force_list[i][1] += random_force(T, dt)
        force_list[i][2] += random_force(T, dt)

def Gauss_random():
    a = random.random()
    b = random.random()
    ksi = np.sqrt(-1*np.log(a))*np.cos(2*np.pi*b) / np.sqrt(5)
    theta = np.sqrt(-1*np.log(b))*np.sin(2*np.pi*a) / np.sqrt(5)

    return ksi, theta

def A(force_list, v_list, dt, k, ksi, theta, i):
    A = (0.5*dt**2 * (force_list[i][k] - gamma*v_list[i][k]) + s*dt**(3/2)*(0.5*ksi + theta/(2*np.sqrt(3))))
    return A


def modify_force_list(f_list, v_list, dt):
    n = len(f_list)
    for i in range(n):
        for k in range(3):
            f_list[i][k] += random_force(dt) - gamma*v_list[i][k]

    return f_list

def random_force(dt):
    a = random.random()
    b = random.random()
    f = s*np.sqrt(-2*np.log(a))*np.cos(2*np.pi*b) / np.sqrt(dt)
    return f

def calc_next_moment(TOTAL_FORCES, part_list, v_list, l, dt, r_cut, moment, kinetic_list, potent_list, vir_list, index):
    n = len(part_list)
    if index:
        f_list = modify_force_list(TOTAL_FORCES[moment], v_list, dt)
    else:
        f_list = TOTAL_FORCES[moment]

    for i in range(n):
        for k in range(3):
            part_list[i][k] = part_list[i][k] + v_list[i][k]*dt + (f_list[i][k])*dt**2 / 2

            if part_list[i][k] > l:
                part_list[i][k] = part_list[i][k] - l
            if part_list[i][k] < 0:
                part_list[i][k] = l + part_list[i][k]


    kinetic, potent, f_list_next, vir = calc_force_energy(part_list, v_list, r_cut, dt, l)
    kinetic_list.append(kinetic)
    potent_list.append(potent)
    vir_list.append(vir)
    if index:
        f_list_next = modify_force_list(f_list_next, v_list, dt)
    TOTAL_FORCES.append(f_list_next)

    for i in range(n):
        v_list[i] = v_list[i] + (f_list[i] + f_list_next[i]) * dt / 2

    return part_list, v_list, TOTAL_FORCES, kinetic_list, potent_list, vir_list

def calc_pressure(n, T, l, vir):
    p = n*T + vir / (3*l**3)
    return p

def print_cords_in_file(part_list, v_list, moment, l):
    moment_file = open('C:/Users/Home/PycharmProjects/project2/data/isoterm4/test_3'
                       '/xyz cor/moment{}.xyz'.format(moment), 'w')
    print('{}'.format(len(part_list)), file=moment_file)
    print('Lattice="{} 0.0 0.0 0.0 {} 0.0 0.0 0.0 {}" Properties=pos:R:3:velo:R:3 Time={}'.format(l, l, l,moment),
          file=moment_file)

    for i in range(len(part_list)):
        state_part = [part_list[i][0], part_list[i][1], part_list[i][2], v_list[i][0], v_list[i][1], v_list[i][2]]
        print('{}'.format(state_part[0]), '{}'.format(state_part[1]), '{}'.format(state_part[2]),
              '{}'.format(state_part[3]), '{}'.format(state_part[4]), '{}'.format(state_part[5]), file=moment_file)

def print_data_in_file(pressure_list, temprature_list, energy_list):
    temp_file = open('C:/Users/Home/PycharmProjects/project2/data/isoterm6/test_1/data test/temprature.txt', 'w')
    pressure_file = open('C:/Users/Home/PycharmProjects/project2/data/isoterm6/test_1/data test/pressure.txt', 'w')
    energy_file = open('C:/Users/Home/PycharmProjects/project2/data/isoterm6/test_1/energy.txt', 'w')
    n = len(pressure_list)
    for i in range(n):
        print(pressure_list[i], file=pressure_file)
        print(temprature_list[i], file=temp_file)
        print(energy_list[i], file=energy_file)



if __name__ == "__main__":
    pass