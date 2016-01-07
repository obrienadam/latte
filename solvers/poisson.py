from equation.fv_equation import FvEquation, DiffusionTerm, AdvectionTerm
import numpy as np

class Poisson(object):
    def __init__(self, grid, **kwargs):
        self.grid = grid

        self.rho = kwargs.get('rho', 1.)
        self.gamma = kwargs.get('gamma', 1.)
        u = kwargs.get('u', np.ones(grid.core_shape))
        v = kwargs.get('v', np.ones(grid.core_shape))
        self.phi, = self.grid.add_cell_fields(kwargs.get('field_name', 'phi'))
        self.bcs = kwargs.get('bcs', {'type': ['f']*4,
                                      'value': range(4)})

        self.adv_term = AdvectionTerm(grid, u, v, self.rho, scheme='upwind')
        self.diff_term = DiffusionTerm(grid, self.gamma)
        self.phi_eqn = FvEquation(self.phi, bcs=self.bcs)

    def solve(self, **kwargs):
        grid = self.grid
        self.phi_eqn == self.adv_term - self.diff_term
        self.phi_eqn.solve()

        return self.grid

if __name__ == '__main__':
    from grid.finite_volume import FvEquidistantGrid
    import grid.viewers

    g = FvEquidistantGrid(500, 1)

    input = {
        'time_accurate': False,
        'dt': 0.05,
        'max_iters': 10,
        'field_name': 'theta',
        'gamma': 0.1,
        'rho': 1.,
        'u': 2.,
        'v': -1.,
        'bcs': {'type': ['f', 'ng', 'f', 'ng'], 'value': [1, 0, 1, 0]},
    }

    solver = Poisson(g, **input)
    solution = solver.solve()

    grid.viewers.display_fv_solution(solution, input['field_name'], show=True, mark_cell_centers=False, show_grid=False)

