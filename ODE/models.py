import numpy as np
import matplotlib.pyplot as plt

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
        return 0.5 * self.rho_g * np.pi**2 * self.a**2 * (
            12 * self.mu_g / self.rho_g / self.a * u
            + 6 * u**2 / (1 + np.sqrt(2 * self.rho_g * u * self.a / self.mu_g))
            + 0.4 * u**2
        )
    
    def f(self, u):
        """ The function of the rhs in the ODE equation """
        return self.g - self.D(u) / self.m_p
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=3, figsize=(10, 10), sharex=True)
        axes[0].plot(time, u, 'k')
        axes[0].set_xlabel('$t$ [s]')
        axes[0].set_ylabel('$v$ [m/s]')
        axes[1].plot(time, self.Re(u), 'k')
        axes[1].set_ylabel('$Re$')
        axes[2].plot(time, self.C_D(u), 'k')
        axes[2].set_ylabel('$C_D$')
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')

class RHSSquare:
    def f(self, u):
        return - u**2
    
    def u_exact(self, time):
        return 1 / (1 + time)
    
    def ax_plot(self, axes, time, u, u_exact, linestyle):
        axes[0].plot(time, u, linestyle)
        axes[0].set_xlabel('$t$ [s]')
        axes[1].plot(time, np.abs(u - u_exact), linestyle)
    
    def plot(self, time, u, figtitle, figname):
        fig, axes = plt.subplots(nrows=2, sharex=True)
        u_exact = self.u_exact(time)
        axes[0].plot(time, u_exact, 'k--')
        self.ax_plot(axes, time, u, u_exact, 'k')
        axes[0].legend(['Exact', 'Simulation'])
        fig.suptitle(figtitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.97])
        fig.savefig(figname, bbox_inches='tight')