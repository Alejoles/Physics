# p5.py

################################################################################
# LIBRERIAS IMPORTADAS Y VARIABLES GLOBALES

import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GM = 4 * (np.pi)**2

################################################################################
# DEFINICION DE CLASES Y FUNCIONES

def universal_gravity(state, params):
    x, y, vx, vy, _ = state
    mass, gu = params
    r2 = x**2 + y**2
    r3 = np.sqrt(r2) * r2
    ax = -gu * x / r3
    ay = -gu * y / r3
    return vx, vy, ax, ay, 1.

def area(deltat, x, y, vx, vy):
    a = (1/2) * deltat * (x * vy - y * vx)
    return a


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES
k = 2 * np.pi
x0, y0, v0, a0, m = 1., .0, k, 90., 1.
sim_params = m, GM
deltat = 0.01/256

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

earth = pt.Particle("Earth", x0, y0, v0, a0, m)
print(earth)
earth_force = fr.Forces(universal_gravity, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, m2, deltat)
xpos, ypos, tpos = [], [], []
vpos = []
cont = 0

areas = []
for i in range(200000):
    xc, yc, vxc, vyc, tc = earth.get_state()
    a = round(area(deltat, xc, yc, vxc, vyc), 5)
    areas.append(a)
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    euler.do_step()
    xc1, yc1, vxc1, vyc1, tc1 = earth.get_state()
    if yc < 0 and yc1 >= 0:
        cont += 1

print("Orbitas:", cont)
print("Area Ini:", areas[0], ", Area Fin:", areas[-1])
# print(areas)
# print(tpos)

# fig, ax = plt.subplots()
# ax.plot(xpos, ypos, '--', label=m2)
#
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Orbits:' + str(cont))



fig, ax = plt.subplots()
ax.plot(tpos, areas, '--', label=m2)

ax.set(xlabel='t (yr)', ylabel="A (a.u.**2)",
       title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Area Circ')

ax.grid()
fig = plt.gcf()
plt.legend()
fig.set_size_inches(10., 5.)
plt.savefig("p5areaCirc.jpg", dpi = 100)
plt.show()


# Para orbita circular x0 = 1, y0 = 0, vx = 0, vy = 2pi, deltat = 0.01/256
# Se cumple la segunda ley de kepler

# Para orbita eliptica x0 = 1.5, y0 = 0, vx = 0, vy = 1.5, deltat = 0.01/256
# Se cumple la segunda ley de kepler
