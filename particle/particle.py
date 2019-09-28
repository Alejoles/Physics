# particle.py

################################################################################
# LIBRERIAS IMPORTADAS Y VARIABLES GLOBALES

import numpy as np

# GRAV = 1.
# DRAG = 4.

GRAV = 1.
DRAG = 0.001 

################################################################################
# DEFINICION DE CLASES Y FUNCIONES

class Particle:

    def __init__(self, label, x0, y0, v0, a0, m0 = 1., t0 = 0.):
        self.label = label
        self.x = x0
        self.y = y0
        self.vx = v0 * np.cos(np.radians(a0))
        self.vy = v0 * np.sin(np.radians(a0))
        self.m = m0
        self.t = t0
        self.force = None


    def __str__(self):
        msg = "Particle: {} \nmass = {:.4f}, t = {:.4f} \n".format(self.label,
                                                                   self.m,
                                                                   self.t)
        msg += "r = ({:.4f}, {:.4f}) \nv = ({:.4f}, {:.4f})".format(self.x,
                                                                    self.y,
                                                                    self.vx,
                                                                    self.vy)
        return msg


    def get_state(self):
        return self.x, self.y, self.vx, self.vy, self.t


    # def euler_step(self, dt):
    #     """ Euler algorithm """
    #     state = self.x, self.y, self.vx, self.vy, self.t
    #     dxdt, dydt, dvxdt, dvydt, dtdt = self.force.get_force(state)
    #     self.x = self.x + dxdt * dt
    #     self.y = self.y + dydt * dt
    #     self.vx = self.vx + dvxdt * dt
    #     self.vy = self.vy + dvydt * dt
    #     self.t = self.t + dtdt * dt


    def euler_step(self, dt):
        """ Euler algorithm """
        self.x = self.x + self.vx * dt
        self.y = self.y + self.vy * dt
        self.vx = self.vx - (DRAG / self.m) * self.vx * dt
        self.vy = self.vy - (DRAG / self.m) * self.vy * dt - GRAV * dt
        self.t = self.t + dt


    def euler_cromer_step(self, dt):
        """ Euler-Cromer algorithm """
        self.vx = self.vx - (DRAG / self.m) * self.vx * dt
        self.vy = self.vy - (DRAG / self.m) * self.vy * dt - GRAV * dt
        self.x = self.x + self.vx * dt
        self.y = self.y + self.vy * dt
        self.t = self.t + dt


    def midpoint_step(self, dt):
        """ Midpoint algorithm """
        vx, vy = self.vx, self.vy
        self.vx = self.vx - (DRAG / self.m) * self.vx * dt
        self.vy = self.vy - (DRAG / self.m) * self.vy * dt - GRAV * dt
        self.x = self.x + .5 * (self.vx + vx) * dt
        self.y = self.y + .5 * (self.vy + vy) * dt
        self.t = self.t + dt


    def set_force(self, netforce):
        self.force = netforce


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

if __name__ == "__main__":
    mars = Particle("Mars", 1., 2., 3., 4.)
    print(mars)
