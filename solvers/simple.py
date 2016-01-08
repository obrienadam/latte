from solver import Solver
from equation.fv_equation import FvEquation, DiffusionTerm, AdvectionTerm, Source

class Simple(Solver):
    def __init__(self, grid, **kwargs):
        super(Simple, self).__init__(grid)

        # Configuration
        self.adv_scheme = kwargs.get('advection_scheme', 'upwind')
        self.omega_momentum = kwargs.get('omega_momentum', 0.7)
        self.omega_p_corr = kwargs.get('omega_p_corr', 0.3)
        self.rho = kwargs.get('rho', 1.)
        self.mu = kwargs.get('mu', 0.01)

        # Initialize fields
        self.u, self.v, self.p, self.p_corr = grid.add_cell_fields('u', 'v', 'p', 'p_corr')
        self.dp_x, self.dp_y = grid.add_cell_fields('dp_x', 'dp_y')
        self.d, = grid.add_cell_fields('d')

        self.uf, self.vf = grid.add_link_fields('uf', 'vf')

        self.u_eqn = FvEquation(self.u)
        self.v_eqn = FvEquation(self.v)
        self.p_corr_eqn = FvEquation(self.p_corr)

        self.laplacian = DiffusionTerm(self.grid, self.mu)

        self._setup_bcs(kwargs.get('bcs', {'type': ['outlet']*4, 'value': [0.]*4}))
        self._compute_interpolation_coeffs(**kwargs)

    def solve(self):
        u_source = Source(self.grid.core_shape)
        v_source = Source(self.grid.core_shape)

        pf = self._interpolate_faces(self.p)
        self.dp_x[1:-1, 1:-1], self.dp_y[1:-1, 1:-1] = self._compute_gradient(pf)

        u_source.b[:, :] = -self.areas*self.dp_x[1:-1, 1:-1] - 1
        v_source.b[:, :] = -self.areas*self.dp_y[1:-1, 1:-1]

        self.u_eqn == AdvectionTerm(self.grid, self.uf, self.vf, self.rho, scheme=self.adv_scheme) - self.laplacian == u_source
        self.v_eqn == AdvectionTerm(self.grid, self.uf, self.vf, self.rho, scheme=self.adv_scheme) - self.laplacian == v_source

        self.u_eqn.relax(self.omega_momentum)
        self.v_eqn.relax(self.omega_momentum)

        self.u_eqn.solve()
        self.v_eqn.solve()

        self.uf[:, :], self.vf[:, :], df = self._rhie_chow_interpolate_faces(self.u, self.v, self.p)

        self.p_corr_eqn == DiffusionTerm(self.grid, self.rho, df)

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

    def _compute_d(self):
        d, = self.grid.get_cell_fields('d')

        d[1:-1, 1:-1] = self.grid.core_cell_areas/0.5*(self.u_eqn.a_core[:, :, 4] + self.v_eqn.a_core[:, :, 4])
        d[0, :] = d[1, :]
        d[-1, :] = d[-2, :]
        d[:, 0] = d[:, 1]
        d[:, -1] = d[:, -2]

        return d

    def _rhie_chow_interpolate_faces(self, u, v, p):
        uf, vf, pf= self._interpolate_faces(u), self._interpolate_faces(v), self._interpolate_faces(p)
        dp_x, dp_y = self.grid.get_cell_fields('dp_x', 'dp_y')

        self.dp_x[1:-1, 1:-1], self.dp_y[1:-1, 1:-1] = self._compute_gradient(pf)

        self.d = self._compute_d()
        df = self._interpolate_faces(self.d)

        sn, rc = self.sn, self.rc
        p = self.p

        den = sn[:, :, :, 0]*rc[:, :, :, 0] + sn[:, :, :, 1]*rc[:, :, :, 1]

        uf[:, :, 0] = uf[:, :, 0] - df[:, :, 0]*(p[2:, 1:-1] - p[1:-1, 1:-1])*sn[:, :, 0, 0]/den[:, :, 0]
        uf[:, :, 1] = uf[:, :, 1] - df[:, :, 1]*(p[1:-1, 2:] - p[1:-1, 1:-1])*sn[:, :, 1, 0]/den[:, :, 1]
        uf[:, :, 2] = uf[:, :, 2] - df[:, :, 2]*(p[:-2, 1:-1] - p[1:-1, 1:-1])*sn[:, :, 2, 0]/den[:, :, 2]
        uf[:, :, 3] = uf[:, :, 3] - df[:, :, 3]*(p[1:-1, :-2] - p[1:-1, 1:-1])*sn[:, :, 3, 0]/den[:, :, 3]

        vf[:, :, 0] = vf[:, :, 0] - df[:, :, 0]*(p[2:, 1:-1] - p[1:-1, 1:-1])*sn[:, :, 0, 1]/den[:, :, 0]
        vf[:, :, 1] = vf[:, :, 1] - df[:, :, 1]*(p[1:-1, 2:] - p[1:-1, 1:-1])*sn[:, :, 1, 1]/den[:, :, 1]
        vf[:, :, 2] = vf[:, :, 2] - df[:, :, 2]*(p[:-2, 1:-1] - p[1:-1, 1:-1])*sn[:, :, 2, 1]/den[:, :, 2]
        vf[:, :, 3] = vf[:, :, 3] - df[:, :, 3]*(p[1:-1, :-2] - p[1:-1, 1:-1])*sn[:, :, 3, 1]/den[:, :, 3]

        return uf, vf, df

if __name__ == '__main__':
    from grid.finite_volume import FvEquidistantGrid
    from grid.viewers import display_fv_solution
    import numpy as np

    bcs = {
        'type': ['outlet', 'wall', 'inlet', 'wall'],
        'value': [0., 0., 1., 0.],
    }

    input = {
        'bcs': bcs,
        'rho': 1.2,
        'mu': 5e-3,
        'omega_momentum': 0.5,
        'omega_p_corr': 0.2,
        'advection_scheme': 'upwind',
    }

    g = FvEquidistantGrid(70, 1)
    simple = Simple(g, **input)

    for i in xrange(400):
        simple.solve()

    u, v = g.get_cell_fields('u', 'v')
    vel, = g.add_cell_fields('vel')

    vel[:, :] = np.sqrt(u*u + v*v)
    print np.max(np.max(vel))
    display_fv_solution(g, 'u', show=True, show_grid=False, mark_cell_centers=False)