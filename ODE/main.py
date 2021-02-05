import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

class FreeFall:
    """ Class of evaluation of FreeFall problem """
    def __init__(self, a, rho_p, rho_g, mu_g, g):
        # Sphere properties
        self.a = a
        self.rho_p = rho_p
        self.m_p = 4 / 3 * np.pi * a**3 * rho_p

        # Properties of the gas
        self.rho_g = rho_g
        self.mu_g = mu_g

        # Gravity coefficient
        self.g = g 
    
    def Re(self, u):
        """ Compute the Reynolds number of the sphere in FreeFall """
        return 2 * self.rho_g * u * self.a / self.mu_g
    
    def C_D(self, u):
        """ Compute the drag coefficient of the sphere """
        return 24 / self.Re(u) + 6 / (1 + np.sqrt(self.Re(u))) + 0.4

    def D(self, u):
        """ Compute the drag force of the sphere in free fall """
        return np.where(u > 0, 0.5 * self.rho_g * np.pi**2 * self.a**2 * u**2 * self.C_D(u), 0)
    
    def f(self, u):
        """ The function of the rhs in the ODE equation """
        return self.g - self.D(u) / self.m_p
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)
        axes[0].plot(time, u, 'k')
        axes[0].set_ylabel('$v$ [m/s]')
        axes[1].plot(time, self.Re(u), 'k')
        axes[1].set_ylabel('$Re$')
        axes[2].plot(time, self.C_D(u), 'k')
        axes[2].set_ylabel('$C_D$')
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')

class ODESim:
    def __init__(self, times, model):
        # Time related variables
        self.times = times
        self.ntimes = len(times)
        self.dts = times[1:] - times[:-1]

        # Model and scheme
        self.model = model

        # variable of interest
        self.v = np.zeros_like(times)
    
    def forwardEuler(self):
        """ Use model to apply forward Euler scheme """
        for i in range(1, self.ntimes):
            self.v[i] = self.v[i - 1] + self.dts[i - 1] * self.model.f(self.v[i - 1])
    
    def midpoint(self):
        """ Use model to apply midpoint formula """
        self.v[1] = self.v[0] + self.dts[0] * self.model.f(self.v[0])
        for i in range(2, self.ntimes):
            self.v[i] = self.v[i - 2] + 2 * self.dts[i - 1] * self.model.f(self.v[i - 1])

if __name__ == '__main__':
    fig_dir = Path('figures/')
    fig_dir.mkdir(exist_ok=True)

    # Time initialization
    tmin, tend, ntimes = 0, 25, 101
    times = np.linspace(tmin, tend, ntimes)
    model = FreeFall(0.01, 917, 0.9, 1.69e-5, 9.81)

    # Class init
    sim = ODESim(times, model)

    # Apply scheme and plot the results
    sim.forwardEuler()
    model.plot(sim.times, sim.v, "Forward Euler", fig_dir.name + '/fwd_euler')

    # Apply midpoint and plot the results
    sim.midpoint()
    model.plot(sim.times, sim.v, "Midpoint formula", fig_dir.name + '/midpoint')