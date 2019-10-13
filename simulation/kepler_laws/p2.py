# p2.py

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

def cond_init(x, y, vx):
    v = np.sqrt(GM / np.sqrt(x**2 + y**2))
    vy = np.sqrt(v**2 - vx**2)
    return float(v), vy

def universal_gravity(state, params):
    x, y, vx, vy, _ = state
    mass, gu = params
    r2 = x**2 + y**2
    r3 = np.sqrt(r2) * r2
    ax = -gu * x / r3
    ay = -gu * y / r3
    return vx, vy, ax, ay, 1.

def total_mechanic_energy(v, m):
    emt = (1/2) * m * (v**2) - GM * m
    return emt

################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

v0, vy = cond_init(1, 0, 0)
x0, y0, v0, a0, m = 1., .0, v0, 90., 1.
sim_params = m, GM
deltat = 0.01/499.998 # para euler cromer deltat = 0.01/499.998

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

earth = pt.Particle("Earth", x0, y0, v0, a0, m) # Revisar por que a0 es 90?
earth_force = fr.Forces(universal_gravity, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, m2, deltat) # El metodo mas estable es Euler Cromer
xpos, ypos, tpos = [], [], []
cont = 0

em = []
for i in range(200000):
    xc, yc, vxc, vyc, tc = earth.get_state()
    v = np.sqrt(vxc**2 + vyc**2)
    emt = total_mechanic_energy(v, m)
    em.append(round(emt, 2))
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    euler.do_step()
    xc1, yc1, vxc1, vyc1, tc1 = earth.get_state()
    if yc < 0 and yc1 >= 0:
        cont += 1
print("Method:", m2, ", dt:", deltat)
for j in em:
    print("Energia Mec Tot:", j)
print("Energia Mec Tot In:", em[0], ", Energia Mec Tot In:", em[-1])
print("Orbitas:", cont)
# print(tpos)

# fig, ax = plt.subplots()
# ax.plot(xpos, ypos, '--', label=m1)
#
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Simulation Planet Motion, deltat: ' + str(deltat))
# ax.grid()

fig, ax = plt.subplots()
ax.plot(tpos, em, '--', label=m2)

ax.set(xlabel='t (yr)', ylabel='emt',
       title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
ax.grid()

plt.legend()
plt.savefig("p2emtConsRound2.jpg")
plt.show()




# El metodo que mejor conserva la energia es Euler Cromer, lo cual es un poco
# logico ya que es el mas estable con un deltat fijo

# La energia no se conserva de forma exacta ya que estamos hablando de una
# trayectoria en forma de elipse, viendo el caso particular de la trayectoria
# circular, si tenemos deltat = 0.01/300.000001, 150000 iter, al principio y
# al final, la emt es aprox la misma con 5dec.
# Ademas con deltat = 0.01/499.998, 200000 iter al principio y al final,
# la emt es aprox la misma con 7dec.
