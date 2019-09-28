import sys
sys.path.insert(0, "../")
print(sys.path)

import particle.particle as pt

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def exact_sol(tf, x0, y0, v0, a0):
    v0x = v0 * np.cos(np.radians(a0))
    v0y = v0 * np.sin(np.radians(a0))
    xa = x0 + v0x * tf
    ya = y0 + v0y * tf - .5 * pt.GRAV * tf ** 2
    return xa, ya

###############################################

deltat = 0.01
x0, y0, v0, a0 = 0., 0., 1., 45.

###############################################

ball = pt.Particle("Ball", x0, y0, v0, a0)
print(ball)

# implement Euler with a = -g
# compare to algebraic solution
# x = x0 + v0x * t
# y = y0 + v0y * t +
# Revisar Fotos para implementar tarea

xpos = []
ypos = []
tpos = []

xana = []
yana = []

x1, y1, vx1, vy1, t1 = ball.get_state()
xpos.append(x1)
ypos.append(y1)
tpos.append(t1)

x2, y2 = exact_sol(t1, x0, y0, v0, a0)

while(y1 >= 0):
    ball.step(deltat)

    x1, y1, vx1, vy1, t1 = ball.get_state()

    xpos.append(x1)
    ypos.append(y1)
    tpos.append(t1)

    # print("Estado de ball en el tiempo t = {:.4f}".format(t1))
    # print("Pos = ({:.4f}, {:.4f})".format(x1, y1))
    # print("Vel = ({:.4f}, {:.4f})".format(vx1, vy1))
    # print()

for t in tpos:
    x2, y2 = exact_sol(t, x0, y0, v0, a0)
    xana.append(x2)
    yana.append(y2)

##################################################

fig, ax = plt.subplots()
# ax.plot(tpos, ypos, "--")
# ax.plot(tpos, yana, "-")
ax.set(xlabel = "Time(t)", ylabel = "ypos", title = "Movement")

errory = []
for l in range(0, len(ypos)):
    yamye = abs(yana[l] - ypos[l])
    errory.append(yamye)

errorx = []
for k in range(0, len(xpos)):
    xamxe = abs(xana[k] - xpos[k])
    errorx.append(xamxe)

ax.plot(tpos, errory, "--")
ax.plot(tpos, errorx, "-")

ax.grid()

plt.show()

# tarea: Volver a integrar con F = mg - cv (como vectores)
# step -> euler_step
# euler_cromer_step cambiar orden, primero velocidad y luego distancia
