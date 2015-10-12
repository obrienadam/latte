class Solver(object):
    def __init__(self, grid):
        self.grid = grid

class Poisson(Solver):
    def __init__(self, grid):
        super(Poisson, self).__init__(grid)
        self.grid.add_fields('phi', 'mu', 'grad_phi_x', 'grad_phi_y')

    def solve(self):
        pass
