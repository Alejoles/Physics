class Solver:

    def __init__(self, objects, nmethod, delta):
        self.step = delta
        self.objs = objects if isinstance(objects, list) else [objects]
        self.algo = Solver.methods[nmethod]

    def __str__(self):
        strng = "Solver state\n"
        strng += "Method = {}".format(self.get_method())
        strng += ", step = {}\n\n".format(self.get_step())
        for idx in range(len(self.objs)):
            strng += "Object[{}]: ".format(idx)
            strng += "{}\n".format(self.objs[idx])
        return strng

    def set_method(self, nmethod):
        self.algo = Solver.methods[nmethod]

    def get_method(self):
        return [n for n, m in Solver.methods.items() if m == self.algo][0]

    def set_step(self, ts):
        self.step = ts

    def get_step(self):
        return self.step

    def do_step(self):
        self.algo(self, self.step)

    def euler_step(self, dt):
        # Euler algorithm
        final_state = []

        for idx in range(len(self.objs)):
            current = self.objs[idx]
            state = current.get_state()

            x, y, vx, vy, t = state
            dxdt, dydt, dvxdt, dvydt, dtdt = current.force.get_force(state)

            x = x + dxdt * dt
            y = y + dydt * dt
            vx = vx + dvxdt * dt
            vy = vy + dvydt * dt
            t = t + dtdt * dt

            final_state.append((x, y, vx, vy, t))

        for idx in range(len(self.objs)):
            self.objs[idx].set_state(*final_state[idx])

    def euler_cromer_step(self, dt):
        # Euler-Cromer algorithm
        state = self.objs.get_state()

        x, y, vx, vy, t = state
        _, _, dvxdt, dvydt, dtdt = self.objs.force.get_force(state)

        vx = vx + dvxdt * dt
        vy = vy + dvydt * dt
        x = x + vx * dt
        y = y + vy * dt
        t = t + dtdt * dt

        self.objs.set_state(x, y, vx, vy, t)

    def midpoint_step(self, dt):
        # Midpoint algorithm
        state = self.objs.get_state()

        x, y, vx, vy, t = state
        dxdt, dydt, dvxdt, dvydt, dtdt = self.objs.force.get_force(state)

        vx = vx + dvxdt * dt
        vy = vy + dvydt * dt
        x = x + .5 * (vx + dxdt) * dt
        y = y + .5 * (vy + dydt) * dt
        t = t + dtdt * dt

        self.objs.set_state(x, y, vx, vy, t)

    methods = {"Euler": euler_step,
               "Euler-Cromer": euler_cromer_step,
               "Midpoint": midpoint_step}


if __name__ == "__main__":
	import sys
	sys.path.insert(0, '../')
    
	import particle.particle as pt

	ball = pt.Particle("aa", 3., 2., 1., 45.)
	euler = Solver(ball, "Midpoint", 0.25)
	print("Integrator", euler)

	euler.set_step(0.01)
	euler.set_method("Euler-Cromer")
	print("Method:", euler.get_method() + "; Step:", euler.get_step())


"""
	ball2 = pt.Particle("aa2", 3., 2., 1., 45.)
	euler2 = Solver(ball2, "Euler-Cromer", 0.25)
	print("Integrator", euler2)

	euler2.set_step(0.01)
	euler2.set_method("Euler")
	print("Method:", euler2.get_method() + "; Step:", euler2.get_step())


"""

"""

	SQRT(GM/R) = V    ||    SQRT[ (4*PI^(2)*R^(3)) / GM] = T
	
	R^2 = X^2 + Y^2

"""








