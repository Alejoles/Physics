# B.py

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

W = 0. # 3% de la ace, asi a = v**2 / r, a = (vx**2 + vy**2) / np.sqrt(x**2 + y**2)
       # a = 20, 3% de 20 es 0.6
W2 = 0.6
print("Con W =", W)
print("Luego de una orbita, W =", W2)


################################################################################
# DEFINICION DE CLASES Y FUNCIONES

def cond_init(x, y, vx):
    v = np.sqrt(GMe / np.sqrt(x**2 + y**2))
    vy = np.sqrt(v**2 - vx**2)
    return float(v), vy

def universal_gravity(state, params):
    x, y, vx, vy, _ = state
    mass, gu, W = params
    r2 = x**2 + y**2
    r3 = np.sqrt(r2) * r2
    ax = (-gu * x / r3) + W
    ay = (-gu * y / r3)
    return vx, vy, ax, ay, 1.

def total_mechanic_energy(v, m, x, y): # cambiar con A & W, con A no cambia,
                                       # con W se agrega -W*x
    emt = (1/2) * m * (v**2) - (GMe * m / np.sqrt(x**2 + y**2)) - W * x   # Julian nos la dio asi, hay que averiguar por que
    return emt


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

v0, vy = cond_init(1, 0, 0)
# print(v0)
x0, y0, v0, a0, m = 1., .0, v0, 90., 1.
sim_params = m, GMe, W
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
                      # y sin W, pero luego se cambia W
            W = W2
            sim_params = m, GMe, W
            satellite_force = fr.Forces(universal_gravity, sim_params)
            satellite.set_force(satellite_force)


print("Method:", m2, ", dt:", deltat)
# for j in em:
#     print("Energia Mec Tot:", j)
print("Energia Mec Tot In:", em[0], ", Energia Mec Tot Fin:", em[-1])
print("Orbitas:", cont)


# emmt = tpos[em.index(min(em))]
# # print(emmt)
# emmx = xpos[tpos.index(emmt)]
# emmy = ypos[tpos.index(emmt)]
# print("Punto donde emt es min:, x = ", round(emmx, 5), "& y =", round(emmy, 5))
# print("Con t =", round(emmt, 4), ", emtm =", min(em))


# GRAFICA xVSy
# fig, ax = plt.subplots()
# ax.plot(xpos, ypos, '--', label=m2)
#
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Simulation Planet Motion with Solar Wind, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
# ax.grid()

# GRAFICA tVSemt
fig, ax = plt.subplots()
ax.plot(tpos, em, '--', label=m2)

ax.set(xlabel='t (h)', ylabel='emt',
       title='Simulation Planet Motion with Solar Wind, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
ax.grid()

# GRAFICA tVSvel
# fig, ax = plt.subplots()
# ax.plot(tpos, vel, '--', label=m2)
#
# ax.set(xlabel='t (h)', ylabel='v(E.U / h)',
#        title='Simulation Planet Motion with Stellar Drag, deltat: ' + str(deltat) + ', Orbits:' + str(cont))
# ax.grid()

plt.legend()
plt.savefig("BemtConW.PNG")
plt.show()

# Igual que con el arrastre estelar, si W != 0, el satelite tiende a salir
# de su orbita original, ya que esta nueva fuerza va en la direccion x, asi
# en un tiempo determinado se produce un choque con la tierra, y aqui la
# simulacion se totea

# Observe que con x0 = 1, y0 = 0, vx = 0, vy = sqrt(20), deltat = 0.01/128, W = 0
# tenemos una orbita circular, grafica orbitaCircSinW

# Ahora con W = 0.6, el satelite hace 7 orbitas, luego se estrella con la
# tierra, grafica orbitaCircConW, ademas de la misma grafica con zoom

# Como la orbita sin W es idealmente circular, la energia mecanica se conserva
# casi de forma exacta, grafica emtSinW

# Con W = 0.6, la energia mecanica cambia en ciertos puntos, grafica emtConW
