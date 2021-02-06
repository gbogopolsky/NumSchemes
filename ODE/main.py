import os
import numpy as np
import matplotlib.pyplot as plt

def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

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

class ODESim:
    def __init__(self, times, schemes, model, init_value, fig_dir=None):
        # Time related variables
        self.times = times
        self.ntimes = len(times)
        self.dts = times[1:] - times[:-1]

        # Model and schemes
        self.model = model
        self.schemes = schemes
        self.nschemes = len(schemes)

        # variable of interest
        self.v = np.zeros((self.nschemes, self.ntimes))
        self.v0 = init_value

        # Figures directory
        if not fig_dir is None:
            self.fig_dir = f'figures/{fig_dir}/'
            create_dir(self.fig_dir)
    
    def forwardEuler(self, v):
        """ Use model to apply forward Euler scheme """
        v[0] = self.v0
        for i in range(1, self.ntimes):
            v[i] = v[i - 1] + self.dts[i - 1] * self.model.f(v[i - 1])
    
    def midpoint(self, v):
        """ Use model to apply midpoint formula """
        v[0] = self.v0
        v[1] = v[0] + self.dts[0] * self.model.f(v[0])
        for i in range(2, self.ntimes):
            v[i] = v[i - 2] + 2 * self.dts[i - 1] * self.model.f(v[i - 1])
    
    def multi_step2(self, v):
        """ Most accurate explicit 2multistep method """
        v[0] = self.v0
        v[1] = v[0] + self.dts[0] * self.model.f(v[0])
        for i in range(2, self.ntimes):
            v[i] = - 4 * v[i - 1] + 5 * v[i - 1] + self.dts[i - 1] * \
                        (4 * self.model.f(v[i - 1]) + 2 * self.model.f(v[i - 2]))
    
    def run_schemes(self):
        """ Apply scheme and plot the results """
        for i_scheme, name_scheme in enumerate(self.schemes):
            scheme = getattr(self, name_scheme)
            scheme(self.v[i_scheme, :])
    
    def plot(self):
        for i_scheme, name_scheme in enumerate(self.schemes):
            self.model.plot(self.times, self.v[i_scheme, :], name_scheme, self.fig_dir + name_scheme)

class ConvergenceODE:
    def __init__(self, tmin, tend, list_nts, schemes, model, init_value, fig_dir):
        self.tmin = tmin
        self.tend = tend
        self.list_nts = list_nts
        self.times = [np.linspace(tmin, tend, ntimes) for ntimes in list_nts]
        self.sims = [ODESim(times, schemes, model, init_value) for times in self.times]
        self.model = model
        self.schemes = schemes
        self.linestyles = ['k-.', 'k--', 'k:']
        self.fig_dir = f'figures/{fig_dir}/'
        create_dir(self.fig_dir)

    def run_convergence(self):
        for sim in self.sims:
            sim.run_schemes()
    
    def plot_errors(self):
        for i_scheme, scheme in enumerate(self.schemes):
            fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
            for i_sim, sim in enumerate(self.sims):
                u_exact = self.model.u_exact(sim.times)
                self.model.ax_plot(axes, sim.times, sim.v[i_scheme, :], u_exact, self.linestyles[i_sim])
            axes[0].plot(sim.times, u_exact, 'k')
            axes[0].legend([rf'$\Delta x$ = {(self.tend - self.tmin) / (nts - 1):.1e}' for nts in self.list_nts] + ['Exact'])
            fig.suptitle(f'Convergence for {scheme}')
            fig.tight_layout(rect=[0, 0.03, 1, 0.97])
            fig.savefig(self.fig_dir + scheme, bbox_inches='tight')


if __name__ == '__main__':
    # Time initialization
    tmin, tend, ntimes = 0, 25, 101
    times = np.linspace(tmin, tend, ntimes)
    model = FreeFall(0.01, 917, 0.9, 1.69e-5, 9.81)
    sim = ODESim(times, ['forwardEuler', 'midpoint', 'multi_step2'], model, 0.0, fig_dir='free_fall/')
    sim.run_schemes()
    sim.plot()

    # Second test with other model
    tmin, tend, ntimes = 0, 10, 101
    times = np.linspace(tmin, tend, ntimes)
    sim = ODESim(times, ['forwardEuler', 'midpoint', 'multi_step2'], RHSSquare(), 1.0, fig_dir='rhs_square')
    sim.run_schemes()
    sim.plot()

    # Convergence test
    tmin, tend = 0, 10
    list_nts = [101, 201, 401]
    cvg_sim = ConvergenceODE(tmin, tend, list_nts, ['forwardEuler', 'midpoint', 'multi_step2'], RHSSquare(), 1.0, 'cvg/')
    cvg_sim.run_convergence()
    cvg_sim.plot_errors()