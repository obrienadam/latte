from grid.finite_volume import FvEquidistantGrid

bcs = {
    'type': ['wall', 'wall', 'wall', 'wall'],
    'value': [0., 1., 0., 0.],
}

input = {
    'bcs': bcs,
    'rho': 1.4,
    'mu': 8.94e-4,
    'omega_momentum': 0.7,
    'omega_p_corr': 0.5,
    'advection_scheme': 'upwind',
    'maxiters': 400,
}

grid = FvEquidistantGrid(120, 1.)

if __name__ == '__main__':
    import solvers.piso as ps
    import grid.viewers as viewers
    import numpy as np

    solver = ps.Piso(grid, **input)
    out_format = 'Solution completion: {}%'

    for iter_no in xrange(input['maxiters']):
        solver.solve()
        print out_format.format((iter_no + 1.)/input['maxiters']*100.)

    vel, = grid.add_cell_fields('vel')

    u, v = grid.get_cell_fields('u', 'v')
    vel[:, :] = np.sqrt(u**2 + v**2)

    viewers.display_fv_solution(grid, 'vel', show=True)