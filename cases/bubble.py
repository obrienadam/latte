from grid.finite_volume import FvRectilinearGrid

bcs = {
    'type': ['wall', 'outlet', 'wall', 'wall'],
    'value': [0., 0., 0., 0.],
}

input = {
    'bcs': bcs,
    'rho_1': 1.,
    'rho_2': 1.,
    'mu_1': 0.001,
    'mu_2': 0.001,
    'omega_momentum': 0.7,
    'omega_p_corr': 0.3,
    'advection_scheme': 'upwind',
    'maxiters': 200,
    'n_p_corrs': 3,
}

grid = FvRectilinearGrid((60, 180), (1, 3))

if __name__ == '__main__':
    from solvers.multiphase_piso import MultiphasePiso
    import grid.viewers as viewers
    import initial_conditions.circle

    solver = MultiphasePiso(grid, **input)

    f, = grid.get_cell_fields('f')
    coords = grid.cell_nodes
    initial_conditions.circle.init_circle(grid, 'f', (0.5, 1), 0.2, 1.)

    out_format = 'Solution completion: {}%'

    """
    for iter_no in xrange(input['maxiters']):
        solver.solve()
        print out_format.format((iter_no + 1.)/input['maxiters']*100.)
    """

    solver._compute_curvature()

    viewers.quiver(grid, 'df_x', 'df_y', show=True)