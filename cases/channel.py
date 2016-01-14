from grid.finite_volume import FvRectilinearGrid

bcs = {
    'type': ['outlet', 'wall', 'inlet', 'wall'],
    'value': [0., 0., 1., 0.],
}

input = {
    'bcs': bcs,
    'rho': 1.,
    'mu': 0.001,
    'omega_momentum': 0.7,
    'omega_p_corr': 0.3,
    'advection_scheme': 'upwind',
    'maxiters': 100,
}

grid = FvRectilinearGrid((120, 40), (3, 1))

if __name__ == '__main__':
    from solvers.piso import Piso
    import grid.viewers as viewers
    import numpy as np

    solver = Piso(grid, **input)

    out_format = 'Solution completion: {}%'

    for iter_no in xrange(input['maxiters']):
        solver.solve()
        print out_format.format((iter_no + 1.)/input['maxiters']*100.)

    vel, = grid.add_cell_fields('vel')

    u, v = grid.get_cell_fields('u', 'v')
    vel[:, :] = np.sqrt(u**2 + v**2)

    viewers.display_fv_solution(grid, 'vel', show=True)