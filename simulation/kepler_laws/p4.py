# p4.py

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

def dist(z1, z2):
    dis = abs(z2 - z1)
    return dis

def universal_gravity(state, params):
    x, y, vx, vy, _ = state
    mass, gu = params
    r2 = x**2 + y**2
    r3 = np.sqrt(r2) * r2
    ax = -gu * x / r3
    ay = -gu * y / r3
    return vx, vy, ax, ay, 1.

def datos_elipse(xpos, ypos):
    xmin = min(xpos)
    xmax = max(xpos)
    ymin = min(ypos)
    ymax = max(ypos)

    sema = dist(xmin, xmax) / 2
    seme = dist(ymin, ymax) / 2

    exc = np.sqrt(sema**2 - seme**2) / sema

    p = np.sqrt(sema ** 3)

    return sema, seme, exc, p

################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES
m = 1.
sim_params = m, GM
deltat = 0.01/1024

m1 = "Euler"
m2 = "Euler-Cromer"
m3 = "Midpoint"

infoelipses = []

x0 = 1.
vy0 = 1.
for i in range(200):
    print("cond ini, x0 & vy0 =", round(x0, 2))
    v0 = vy0
    x0, y0, v0, a0, m = x0, .0, v0, 90., 1.

    earth = pt.Particle("Earth", x0, y0, v0, a0, m)
    # print(earth)
    earth_force = fr.Forces(universal_gravity, sim_params)
    earth.set_force(earth_force)
    euler = sl.Solver(earth, m2, deltat)

    xpos, ypos, tpos = [], [], []
    cont = 0
    for i in range(200000):
        xc, yc, vxc, vyc, tc = earth.get_state()
        # if tc > 0.415:
        #     break
        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)
        euler.do_step()
        xc1, yc1, vxc1, vyc1, tc1 = earth.get_state()
        if yc < 0 and yc1 >= 0:
            cont += 1

    sema, seme, exc, p = datos_elipse(xpos, ypos)

    print("Sema =", sema, ", Seme =", seme)
    print("Exc =", exc, ", p =", p)
    print()

    infoelipses.append((round(x0, 2), round(vy0, 2), round(sema, 5), round(p, 5), cont))

    x0 += 0.008
    vy0 += 0.008

# print(infoelipses)
# print("Orbitas:", cont)
# print(tpos)

semapos = []
ppos = []
for k in infoelipses:
    semapos.append(k[2])
    ppos.append(k[3])
# print(semapos)
# print(ppos)

tl = []
logsemapos = []
logppos = []
for j in range(len(ppos)):
    tl.append((ppos[j] ** 2) / (semapos[j] ** 3))
    logsemapos.append(np.log(semapos[j]))
    logppos.append(np.log(ppos[j]))
# print(tl)

fig, ax = plt.subplots()
ax.plot(ppos, semapos, '--', label=m2)

ax.set(xlabel='sema (a.u.)', ylabel='p (yr)',
       title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Elipses Distintas: 200')


# fig, ax = plt.subplots()
# ax.plot(logppos, logsemapos, '--', label=m2)
#
# ax.set(xlabel='log(a)', ylabel='log(p)',
#        title='Simulation Planet Motion, deltat: ' + str(deltat) + ', Elipses Distintas: 200')

ax.grid()
plt.legend()
plt.savefig("p4pVSaLOG.jpg")
plt.show()

# Con deltat = 0.01/1024, y0 = 0, vx0 = 0
# 1) x0 = 1.05, vy0 = 1.05
# 2) x0 = 1.1, vy0 = 1.1
# 3) x0 = 1.2, vy0 = 1.2
# 4) x0 = 1.3, vy0 = 1.3
# 5) x0 = 1.4, vy0 = 1.4
# 6) x0 = 1.5, vy0 = 1.5
# 7) x0 = 1.6, vy0 = 1.6
# 8) x0 = 1.7, vy0 = 1.7
# 9) x0 = 1.8, vy0 = 1.8
# 10) x0 = 1.9, vy0 = 1.9
# Hay trayectorias elipticas

# Sema y periodo respectivamente, tomados a mano
# 1) sema = 0.5328, p = 0.3889
# 2) sema = 0.5594, p = 0.4184
# 3) sema = 0.6134, p = 0.4804
# 4) sema = 0.6686, p = 0.5467
# 5) sema = 0.7252, p = 0.6175
# 6) sema = 0.7834, p = 0.6935
# 7) sema = 0.8437, p = 0.7750.
# 8) sema = 0.9063, p = 0.8629
# 9) sema = 0.9717, p = 0.9576
# 10) sema = 1.0403, p = 1.0611
