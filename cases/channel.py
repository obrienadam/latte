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
}

grid = FvRectilinearGrid((120, 40), (3, 1))