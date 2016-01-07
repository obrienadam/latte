from equation.fv_equation import FvEquation, DiffusionTerm, AdvectionTerm

class Simple(object):
    def __init__(self, grid, **kwargs):
        self.grid = grid

        self.u, self.v, self.p, self.p_corr = grid.add_cell_fields('u', 'v', 'p', 'p_corr')
        self.u_f, self.v_f = grid.add_link_fields('u_f', 'v_f')

        self.rho = kwargs.get('rho', 1.)
        self.mu = kwargs.get('mu', 0.1)
        self.adv_scheme = kwargs.get('advection_scheme', 'upwind')

        self.u_eqn = FvEquation(self.u)
        self.v_eqn = FvEquation(self.v)
        self.p_corr_eqn = FvEquation(self.p_corr)

        self.laplacian = DiffusionTerm(self.grid, self.mu)

        self._setup_bcs(kwargs.get('bcs', {'type': ['outlet']*4, 'value': [0.]*4}))
        self._compute_interpolation_coeffs(**kwargs)

    def solve(self):
        self._interpolate_faces(self.u, self.u_f)
        self._interpolate_faces(self.v, self.v_f)

        self.u_eqn == AdvectionTerm(self.grid, self.u_f, self.v_f, self.rho, scheme=self.adv_scheme) - self.laplacian
        self.v_eqn == AdvectionTerm(self.grid, self.u_f, self.v_f, self.rho, scheme=self.adv_scheme) - self.laplacian

        self.u_eqn.solve()
        self.v_eqn.solve()

    def _setup_bcs(self, bcs):
        self.bcs = bcs

        u_bcs = {'type': ['fixed']*4, 'value': [0.]*4}
        v_bcs = {'type': ['fixed']*4, 'value': [0.]*4}
        p_bcs = {'type': ['fixed']*4, 'value': [0.]*4}

        types, values = self.bcs['type'], self.bcs['value']
        u_types, u_values = u_bcs['type'], u_bcs['value']
        v_types, v_values = v_bcs['type'], v_bcs['value']
        p_types, p_values = p_bcs['type'], p_bcs['value']

        for i in xrange(4):
            if types[i] is 'inlet':
                u_types[i] = 'fixed'
                v_types[i] = 'fixed'
                p_types[i] = 'normal_gradient'
                u_values[i], v_values[i] = (values[i], 0.) if i%2 is 0 else (0., values[i])
                p_values[i] = 0.
            elif types[i] is 'outlet':
                u_types[i] = 'normal_gradient'
                v_types[i] = 'normal_gradient'
                p_types[i] = 'fixed'
                u_values[i], v_values[i] = 0., 0.
                p_values[i] = values[i]
            elif types[i] is 'wall':
                u_types[i] = 'fixed'
                v_types[i] = 'fixed'
                p_types[i] = 'normal_gradient'
                u_values[i], v_values[i] = (0., values[i]) if i%2 is 0 else (values[i], 0.)
            else:
                raise ValueError

        self.u_eqn.bcs = u_bcs
        self.v_eqn.bcs = v_bcs
        self.p_corr_eqn.bcs = p_bcs

    def _compute_interpolation_coeffs(self, **kwargs):
        grid = self.grid
        self.alpha, = grid.add_link_fields('alpha')
        alphas = self.alpha

        if kwargs.get('interpolation_method', 'volume_weighted') == 'volume_weighted':
            areas = grid.core_cell_areas

            alphas[:-1, :, 0] = areas[1:, :]/(areas[1:, :] + areas[:-1, :])
            alphas[:, :-1, 1] = areas[:, 1:]/(areas[:, 1:] + areas[:, :-1])
            alphas[1:, :, 2] = areas[:-1, :]/(areas[:-1, :] + areas[1:, :])
            alphas[:, 1:, 3] = areas[:, :-1]/(areas[:, :-1] + areas[:, 1:])
        else:
            raise ValueError

    def _interpolate_faces(self, cell_field, face_field):
        alphas = self.alpha
        face_field[:, :, 0] = alphas[:, :, 0]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 0])*cell_field[2:, 1:-1]
        face_field[:, :, 1] = alphas[:, :, 1]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 1])*cell_field[1:-1, 2:]
        face_field[:, :, 2] = alphas[:, :, 2]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 2])*cell_field[:-2, 1:-1]
        face_field[:, :, 3] = alphas[:, :, 3]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 3])*cell_field[1:-1, :-2]

if __name__ == '__main__':
    from grid.finite_volume import FvEquidistantGrid
    from grid.viewers import display_fv_solution, plot_line
    import numpy as np

    bcs = {
        'type': ['inlet', 'inlet', 'outlet', 'wall'],
        'value': [1., 1., 0., 0.],
    }

    g = FvEquidistantGrid(50, 1)
    simple = Simple(g, bcs=bcs)

    for i in xrange(5):
        simple.solve()

    u, v = g.get_cell_fields('u', 'v')
    vel, = g.add_cell_fields('vel')

    vel[:, :] = np.sqrt(u*u + v*v)

    plot_line(g, 'u', 25, axis=0, show=True)
    display_fv_solution(g, 'vel', show=True)