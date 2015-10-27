from .solver import Solver

class Simple(Solver):
    def __init__(self, grid, **kwargs):
        super(Simple, self).__init__(grid, **kwargs)

        grid.add_fields('u', 'v', 'p', 'm', 'pCorr', cell_centered=True, face_centered=True)
        grid.add_fields('gradP_x', 'gradP_y', 'gradPCorr_x', 'gradPCorr_y', cell_centered=True)

        self.rho = kwargs.get('rho', 1.)
        self.mu = kwargs.get('mu', 0.1)

    def solve(self, progress_bar):

        for iterNo in xrange(self.max_iters):

            progress_bar.setValue(iterNo + 1)




