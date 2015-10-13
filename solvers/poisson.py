from .solver import Solver
import time

class Poisson(Solver):
    def __init__(self, grid, **kwargs):
        super(Poisson, self).__init__(grid, **kwargs)
        self.grid.add_fields('phi', 'mu', cell_centered=True, face_centered=True)

    def solve(self, progress_bar):
        for i in xrange(1, self.max_iters):
            progress_bar.setValue(i)
            time.sleep(0.1)


