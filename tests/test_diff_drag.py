# test_air_drag.py

################################################################################
# LIBRERIAS IMPORTADAS Y VARIABLES GLOBALES

import sys
sys.path.insert(0, "../")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import particle.particle as pt


################################################################################
# DEFINICION DE CLASES Y FUNCIONES

def analytic_solution(tf, v0, a0):
    """This expression assumes x0 = 0 and y0 = 0 and m = 1"""
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = v0x * (1. - np.exp(-pt.DRAG * tf)) / pt.DRAG
    ya = (v0y + pt.GRAV / pt.DRAG) * (1. - np.exp(-pt.DRAG * tf)) \
       - pt.GRAV * tf
    return xa, ya / pt.DRAG

################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

deltat = 0.005
x0, y0, v0, a0 = 0., 0., 1., 75.

ball = pt.Particle("Ball", x0, y0, v0, a0)
ball1 = pt.Particle("Ball1", x0, y0, v0, a0)
ball2 = pt.Particle("Ball2", x0, y0, v0, a0)
print(ball)
print(ball1)
print(ball2)

xpos = []
ypos = []
tpos = []

while True:
    x, y, _, _, t = ball.get_state()
    if y < 0: break
    xpos.append(x)
    ypos.append(y)
    tpos.append(t)
    ball.midpoint_step(deltat)

xana = []
yana = []

xpos1 = []
ypos1 = []

xpos2 = []
ypos2 = []

for t in tpos:
    pt.DRAG = pt.GRAV/2
    x, y = analytic_solution(t, v0, a0)
    if y < 0: break
    xana.append(x)
    yana.append(y)

    x1, y1, _, _, t1 = ball1.get_state()
    t1 = t
    xpos1.append(x1)
    ypos1.append(y1)
    ball1.midpoint_step(deltat)


for t in tpos:
    pt.DRAG = (pt.GRAV/2) ** 2

    x2, y2, _, _, t2 = ball2.get_state()
    t2 = t

    if y2 < 0: break

    xpos2.append(x2)
    ypos2.append(y2)
    ball2.midpoint_step(deltat)

##############################################

fig, ax = plt.subplots()
# ax.plot(tpos, errorx, '--', label = '|x_ana - x_pos|')
# ax.plot(tpos, errory, '-.', label = '|y_ana - y_pos|')

##############################################

ax.plot(xpos, ypos, '-', label = 'Free Falling')
ax.plot(xpos1, ypos1, '-', label = 'Linear Drag')
ax.plot(xpos2, ypos2, '-', label = 'Quadratic Drag')
ax.plot(xana, yana, '--', label = 'Analytic Linear')

##############################################

# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Projectile motion. Method: euler_step')
# ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
#        title='Projectile motion. Method: euler_cromer_step')
ax.set(xlabel='x (a.u.)', ylabel='y (a.u.)',
       title='Different Drags (dt = 0.005, a0 = 75, C = g/2, midpoint)')

ax.grid()

plt.legend()
plt.savefig("DiffDragsVs.png")
plt.show()
