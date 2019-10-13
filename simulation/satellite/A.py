# A.py

################################################################################
# LIBRERIAS IMPORTADAS Y VARIABLES GLOBALES

import sys
sys.path.insert(0, '../../')

import numpy as np
import matplotlib.pyplot as plt

import particle.particle as pt
import solver.solver as sl
import forces.forces as fr

GMe = 20

A = 0. # con 0 <= A <= 0.3538, de 11 a 7 orbitas, luego se totea por errores numericos
A2 = 0.3538
print("Con A =", A)
print("Luego de una orbita, A =", A2)


################################################################################
# DEFINICION DE CLASES Y FUNCIONES

def cond_init(x, y, vx):
    v = np.sqrt(GMe / np.sqrt(x**2 + y**2))
    vy = np.sqrt(v**2 - vx**2)
    return float(v), vy

def universal_gravity(state, params):
    x, y, vx, vy, _ = state
    mass, gu, A = params
    r2 = x**2 + y**2
    r3 = np.sqrt(r2) * r2
    ax = (-gu * x / r3) + A * np.sqrt(vx**2 + vy**2) * x
    ay = (-gu * y / r3) - A * np.sqrt(vx**2 + vy**2) * y
    return vx, vy, ax, ay, 1.

def total_mechanic_energy(v, m, x, y): # cambiar con A & W
    emt = (1/2) * m * (v**2) - GMe * m / np.sqrt(x**2 + y**2)  # Julian nos la dio asi, hay que averiguar por que
    return emt


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

v0, vy = cond_init(1, 0, 0)
# print(v0)
x0, y0, v0, a0, m = 1., .0, v0, 90., 1.
sim_params = m, GMe, A
deltat = 0.01/128

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

satellite = pt.Particle("Satellite", x0, y0, v0, a0, m) # Revisar por que a0 es 90?
satellite_force = fr.Forces(universal_gravity, sim_params)
satellite.set_force(satellite_force)
euler = sl.Solver(satellite, m2, deltat) # El metodo mas estable es Euler Cromer
xpos, ypos, tpos = [], [], []
cont = 0

em = []
vel = []
for i in range(200000):
    xc, yc, vxc, vyc, tc = satellite.get_state()
    v = np.sqrt(vxc**2 + vyc**2)
    vel.append(v)
    emt = total_mechanic_energy(v, m, xc, yc)
    em.append(round(emt, 2))
    xpos.append(xc)
    ypos.append(yc)
    tpos.append(tc)
    euler.do_step()
    xc1, yc1, vxc1, vyc1, tc1 = satellite.get_state()
    if yc < 0 and yc1 >= 0:
        cont += 1
        if cont == 1: # Observe que se permite que se realice una orbita con las condiciones iniciales
                      # y sin A, pero luego se cambia A
            A = A2
            sim_params = m, GMe, A
            satellite_force = fr.Forces(universal_gravity, sim_params)
            satellite.set_force(satellite_force)

print("Method:", m2, ", dt:", deltat)
# for j in em:
#     print("Energia Mec Tot:", j)
print("Energia Mec Tot In:", em[0], ", Energia Mec Tot Fin:", em[-1])
print("Orbitas:", cont)

# emmt = tpos[em.index(max(em))]
# # print(emmt)
# emmx = xpos[tpos.index(emmt)]
# emmy = ypos[tpos.index(emmt)]
# print("Punto donde emt es max:, x = ", round(emmx, 5), "& y =", round(emmy, 5))
# print("Con t =", round(emmt, 4), ", emtm =", max(em))

# GRAFICA xVSy
# fig, ax = plt.subplots()
# ax.plot(xpos, ypos, '--', label=m2)
#
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Simulation Planet Motion with Stellar Drag, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
# ax.grid()

# GRAFICA tVSemt
fig, ax = plt.subplots()
ax.plot(tpos, em, '--', label=m2)

ax.set(xlabel='t (h)', ylabel='emt',
       title='Simulation Planet Motion with Stellar Drag, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
ax.grid()

# GRAFICA tVSvel
# fig, ax = plt.subplots()
# ax.plot(tpos, vel, '--', label=m2)
#
# ax.set(xlabel='t (h)', ylabel='v(E.U / h)',
#        title='Simulation Planet Motion with Stellar Drag, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
# ax.grid()

plt.legend()
plt.savefig("AtVSemt.PNG")
plt.show()

# Observe que con x0 = 1, y0 = 1, vx = 0, vy = sqrt(20), deltat = 0.01/128, A = 0
# tenemos una orbita circular, grafica orbitaCircSinA

# Ahora con A = 0.3538, el satelite hace 6 orbitas, luego se estrella con la
# tierra, grafica orbitaCircConA

# El metodo mas estable es Euler Cromer, por ende fue elegido para realizar las
# modificaciones con Arrastre Estelar

# Sabemos que la velocidad es tangencial a la orbita del satelite, y que la
# fuerza de arrastre es opuesta a la velocidad, asi esta fuerza tiende a frenar
# el satelite, scandolo de su orbita generando un choque en un tiempo determinado
# con el planeta sobre el cual esta orbitando, en este caso La Tierra.

# Graficas tVSemt, tVSvel

# Para ver como cambian los resultados dependiendo el integrador numerico usado,
# ver graficas: resulEuler, resulMidpoint.

# resultEuler: Para el metodo euler, la simulacion se totea completamente con
# A = 0.3536, asi que hay que bajar A, ahora A = 0.16

# resultMidpoint: Para el metodo Midpoint, simulacion se totea completamente con
# A = 0.3536, asi que hay que bajar A, ahora A = 0.1
