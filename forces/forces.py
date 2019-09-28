# forces.py

################################################################################
# LIBRERIAS IMPORTADAS Y VARIABLES GLOBALES

################################################################################
# DEFINICION DE CLASES Y FUNCIONES

class Forces:

    def __init__(self, user_force, simul_params):
        self.force = user_force
        self.params = simul_params
        self.state = None


    def get_state(self):
        return self.state


    def get_force(self, state):
        self.state = state
        return self.force(state, self.params)


################################################################################
# BLOQUE PRINCIPAL DE INSTRUCCIONES

if __name__ == "__main__":
    G, M = 9.8, 7.

    def falling_particle(state, params):
        # Constants like Gravity, Friction are global
        # args = x, y, vx, vy, t
        # returns : dxdt, dydt, dvxdt, dvydt, dtdt
        xc, yc, vxc, vyc, tc = state
        gg, mm = params
        return vxc, vyc, 0., -gg, 1.

    grav_force = Forces(falling_particle, (G, M))
    print("STATE:", grav_force.get_state())

    xi, yi, vxi, vyi, t = 1., 2., 3., 4., 0.
    params_tuple = xi, yi, vxi, vyi, t
    state = grav_force.get_force(params_tuple)
    print(grav_force.get_force, "\nOUTPUT:", state)
    print("STATE:", grav_force.get_state())
