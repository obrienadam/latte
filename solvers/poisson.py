from .solver import Solver

class Poisson(Solver):
    def __init__(self, grid):
        super(Poisson, self).__init__(grid)
        self.grid.add_fields('phi', 'mu', cell_centered=True, face_centered=True)

    def solve(self):
        pass
