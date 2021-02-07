import numpy as np
import matplotlib.pyplot as plt
from .utils import create_dir

class ODESim:
    def __init__(self, times, schemes, model, init_value, fig_dir=None):
        # Time related variables
        self.times = times
        self.ntimes = len(times)
        self.dts = times[1:] - times[:-1]
        self.dt = self.dts[0]

        # Model and schemes
        self.model = model
        self.schemes = schemes
        self.nschemes = len(schemes)

        # variable of interest
        self.v = np.zeros((self.nschemes, self.ntimes, self.model.nd))
        self.v0 = init_value

        # Figures directory
        if not fig_dir is None:
            self.fig_dir = f'figures/{fig_dir}/'
            create_dir(self.fig_dir)
    
    def forwardEuler(self, v):
        """ Use model to apply forward Euler scheme """
        v[0] = self.v0
        for i in range(1, self.ntimes):
            v[i] = v[i - 1] + self.dts[i - 1] * self.model.f(v[i - 1], self.times[i - 1])
    
    def midpoint(self, v):
        """ Use model to apply midpoint formula """
        v[0] = self.v0
        v[1] = v[0] + self.dts[0] * self.model.f(v[0], self.times[0])
        for i in range(2, self.ntimes):
            v[i] = v[i - 2] + 2 * self.dts[i - 1] * self.model.f(v[i - 1], self.times[i - 1])
    
    def multi_step2(self, v):
        """ Most accurate explicit 2multistep method """
        v[0] = self.v0
        v[1] = v[0] + self.dts[0] * self.model.f(v[0], self.times[0])
        for i in range(2, self.ntimes):
            v[i] = - 4 * v[i - 1] + 5 * v[i - 1] + self.dts[i - 1] * \
                        (4 * self.model.f(v[i - 1], self.times[i - 1]) + 2 * self.model.f(v[i - 2], self.times[i - 2]))
    
    def run_schemes(self):
        """ Apply scheme and plot the results """
        for i_scheme, name_scheme in enumerate(self.schemes):
            scheme = getattr(self, name_scheme)
            scheme(self.v[i_scheme, :])
    
    def plot(self):
        for i_scheme, name_scheme in enumerate(self.schemes):
            self.model.plot(self.times, self.v[i_scheme, :], name_scheme, self.fig_dir + name_scheme)
