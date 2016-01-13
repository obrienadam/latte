from grid.finite_volume import FvEquidistantGrid

bcs = {
    'type': ['wall', 'wall', 'wall', 'wall'],
    'value': [0., 1., 0., 0.],
}

input = {
    'bcs': bcs,
    'rho': 1.,
    'mu': 0.001,
    'omega_momentum': 0.3,
    'omega_p_corr': 0.3,
    'advection_scheme': 'upwind',
}

grid = FvEquidistantGrid(60, 1.)