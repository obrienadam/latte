from .solver import Solver
import numpy as np

class Poisson(Solver):
    def __init__(self, grid, **kwargs):
        super(Poisson, self).__init__(grid, **kwargs)
        self.grid.add_fields('phi', 'mu', cell_centered=True, face_centered=True)

    def solve(self, progress_bar):
        for i in xrange(1, self.max_iters):
            x = np.linspace(0, i, 1000)*np.transpose(np.linspace(0, 2*i, 1000))
            progress_bar.setValue(i)


