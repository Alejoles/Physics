import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../')

import particle.particle as pt
import forces.forces as fr
import solver.solver as sol


GM = 4. * np.pi**2    # gravitational constant times big mass
NITER = 100000        # max number of iterations


def newton_again(state, params):
    xp, yp, vxp, vyp, _ = state
    gm_ = params
    rp3 = np.power(xp**2 + yp**2, 1.5)
    axp, ayp = -gm_ * xp / rp3, -gm_ * yp / rp3
    return vxp, vyp, axp, ayp, 1


def integrate(obj):
    senel, orbit = 0, 0
    xpos, ypos, tpos = [], [], []

    while senel < NITER:
        xc, yc, _, _, tc = obj.get_state()
        if xc > 0 and abs(yc) < 0.001:
            print("TURN", orbit, xc, yc)
            orbit += 1

        xpos.append(xc)
        ypos.append(yc)
        tpos.append(tc)

        obj.do_step()
        senel += 1

    return xpos, ypos, tpos


# BEGINNING-OF-EXECUTION
deltat, num_method = .001, "Euler-Cromer"
x0, y0, v0, a0 = 1., .0, .9, 90
r0 = np.sqrt(x0**2 + y0**2)
v0 = 1.3 * np.sqrt(GM / r0)

# create free-falling Particle
sim_params = GM
planet = pt.Particle("Earth",x0, y0, v0, a0)

the_force = fr.Forces(newton_again, sim_params)
planet.set_force(the_force)

intor = sol.Solver(planet, num_method, deltat)
xvac, yvac, tvac = integrate(intor)

print("x", len(xvac), "y", len(xvac))

# generate plots
fig, ax = plt.subplots()
label = 'x0 = %.2f, y0 = %.2f, v0 = %.2f, a0 = %.2f' % (x0, y0, v0, a0)
ax.plot(xvac, yvac, '-', label=label)

ax.set_aspect('equal')
ax.set(xlabel='x (AU)', ylabel='y (AU)', title='Real free falling')
ax.grid()

plt.legend()
plt.show()
# END-OF-EXECUTION
