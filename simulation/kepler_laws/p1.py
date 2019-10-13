# p1.py

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


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

v0, vy = cond_init(1, 0, 0)
x0, y0, v0, a0, m = 1., .0, v0, 90., 1.
sim_params = m, GM
deltat = 0.01/512

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

earth = pt.Particle("Earth", x0, y0, v0, a0, m)
earth_force = fr.Forces(universal_gravity, sim_params)
earth.set_force(earth_force)
euler = sl.Solver(earth, m1, deltat)
xpos, ypos, tpos = [], [], []
cont = 0



for i in range(200000):
    xc, yc, _, _, tc = earth.get_state()
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    euler.do_step()
    xc1, yc1, _, _, tc1 = earth.get_state()
    if yc < 0 and yc1 >= 0:
        cont += 1

print("Orbitas:", cont) # Aqui nos damos cuenta cuantas orbitas circulares hace el planeta

# print(tpos)
fig, ax = plt.subplots()
ax.plot(xpos, ypos, '--', label=m1)

ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
ax.grid()

plt.legend()
plt.savefig("p1EulerEstable.PNG")
plt.show()



# Revisar por que a0 es 90?

# El metodo mas estable es Euler Cromer

# Para que Euler haga orbitas circulas por muchos periodos, deltat = 0.01/512 aprox
