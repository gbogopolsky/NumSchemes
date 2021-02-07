import numpy as np
import matplotlib.pyplot as plt
from odesolver.solver import ODESim
from odesolver.utils import make_times

class StiffProblem:
    def __init__(self, lambda_1, lambda_2):
        self.nd = 1
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
    
    def f(self, u, t):
        return - self.lambda_1 * u + self.lambda_1 / 10 * np.sin(self.lambda_2 * t)
    
    def plot(self, time, u, figtitle, figname):
        fig, ax = plt.subplots()
        ax.plot(time, u, 'k')
        ax.set_xlabel('$t$ [s]')
        fig.suptitle(figtitle)
        fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    model = StiffProblem(1000, 1)
    tmin, tend = 0, 10
    dts = [1.0e-3, 1.9e-3, 2.0e-3, 2.1e-3]
    for i, dt in enumerate(dts):
        times = make_times(tmin, tend, dt)
        sim = ODESim(times, ['forwardEuler'], model, 1.0, fig_dir=f'stiffproblem/case_{i:d}/')
        sim.run_schemes()
        sim.plot()