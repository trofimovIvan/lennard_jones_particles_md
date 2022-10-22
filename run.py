from test_1 import *

""" Here some code to run simulations"""

l = 9.28
temprature = 0.0
n = 0.1
r_cut = 2.7
sigma = 1




particles_list = []
particles_list = spawn_moleculs(particles_list, l, n*l**3, sigma)
N = len(particles_list)

v_list = []
v_list = get_velocity_to_particles(particles_list, temprature)
v_list = get_norm_molec_v(v_list)

time_list = []
energy_list = []
kinetic_list = []
potent_list = []
temp_list = []
pressure_list = []
vir_list = []

TOTAL_FORCES = []



time_end = 3.0
time = 0.0
dt = 0.001
moment = 0
index = True

while time <= time_end:
    print("time = ", time)
    if moment == 0:
        kinetic, potent, force_list, vir = calc_force_energy(particles_list, v_list, r_cut, dt, l)
        TOTAL_FORCES.append(force_list)
        vir_list.append(vir)

    particles_list, v_list, TOTAL_FORCES, kinetic_list, potent_list, vir_list = calc_next_moment(TOTAL_FORCES, particles_list,
                                              v_list, l, dt, r_cut, moment, kinetic_list, potent_list, vir_list, index)

    temp = kinetic_list[moment]*2/(3*N)

  #  print_cords_in_file(particles_list, v_list, moment, l)

    pressure_list.append(calc_pressure(n, temp, l, vir_list[moment]))

    energy_list.append(kinetic_list[moment] + potent_list[moment])
    time_list.append(time)
    temp_list.append(temp)

    print("temp = ", temp)
    print("energy = ", kinetic_list[moment]+potent_list[moment])
    print('pressure = ', pressure_list[moment])
    moment += 1
    time += dt
    if time >= 2.0:
        index = False


plt.plot(time_list, energy_list, color='blue')
plt.plot(time_list, kinetic_list, color='green')
plt.plot(time_list, potent_list, color='red')
plt.grid()
plt.show()

plt.plot(time_list, temp_list)
plt.show()

plt.errorbar(time_list, pressure_list, fmt='.')
plt.show()

print_data_in_file(pressure_list, temp_list, energy_list)
