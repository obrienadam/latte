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
        self.dp_corr_x, self.dp_corr_y = grid.add_cell_fields('dp_corr_x', 'dp_corr_y')
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

        u_source.b[:, :] = -self.areas*self.dp_x[1:-1, 1:-1]
        v_source.b[:, :] = -self.areas*self.dp_y[1:-1, 1:-1]

        self.u_eqn == (AdvectionTerm(self.grid, self.uf, self.vf, self.rho, scheme=self.adv_scheme) - self.laplacian == u_source)
        self.v_eqn == (AdvectionTerm(self.grid, self.uf, self.vf, self.rho, scheme=self.adv_scheme) - self.laplacian == v_source)

        self.u_eqn.relax(self.omega_momentum)
        self.v_eqn.relax(self.omega_momentum)

        self.u_eqn.iterative_solve(maxiter=20)
        self.v_eqn.iterative_solve(maxiter=20)

        self.uf[:, :], self.vf[:, :], df = self._rhie_chow_interpolate_faces(self.u, self.v, self.p)

        mass_source = Source(self.grid.core_shape)
        mass_source.b[:, :] = self.rho*np.sum((self.uf*self.grid.core_face_norms[:, :, :, 0] + self.vf*self.grid.core_face_norms[:, :, :, 1]), dtype=float, axis=-1)

        self.p_corr_eqn == (DiffusionTerm(self.grid, df, self.rho) == mass_source)
        self.p_corr_eqn.iterative_solve(maxiter=100, tol=1e-4)

        self._correct(df)

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
                p_values[i] = 0.
            elif types[i] is 'wall':
                u_types[i] = 'fixed'
                v_types[i] = 'fixed'
                p_types[i] = 'normal_gradient'
                u_values[i], v_values[i] = (0., values[i]) if i%2 is 0 else (values[i], 0.)
                p_values[i] = 0.
            else:
                raise ValueError

        self.u_eqn.bcs = u_bcs
        self.v_eqn.bcs = v_bcs
        self.p_corr_eqn.bcs = p_bcs

    def _compute_d(self):
        d = self.d

        d[1:-1, 1:-1] = self.grid.core_cell_areas/self.u_eqn.a_core[:, :, 4]
        d[0, :] = d[1, :]
        d[-1, :] = d[-2, :]
        d[:, 0] = d[:, 1]
        d[:, -1] = d[:, -2]

        return d

    def _rhie_chow_interpolate_faces(self, u, v, p):
        pf = self._interpolate_faces(p)
        dp_x, dp_y = self.dp_x, self.dp_y

        dp_x[1:-1, 1:-1], dp_y[1:-1, 1:-1] = self._compute_gradient(pf)

        d = self._compute_d()
        df = self._interpolate_faces(d)

        sn, rc = self.sn, self.rc
        p = self.p

        den = sn[:, :, :, 0]*rc[:, :, :, 0] + sn[:, :, :, 1]*rc[:, :, :, 1]
        p_curr = p[1:-1, 1:-1]

        hu = u + d*dp_x
        hv = v + d*dp_y

        uf, vf = self._interpolate_faces(hu), self._interpolate_faces(hv)

        uf[:, :, 0] = uf[:, :, 0] - df[:, :, 0]*(p[2:, 1:-1] - p_curr)*sn[:, :, 0, 0]/den[:, :, 0]
        uf[:, :, 1] = uf[:, :, 1] - df[:, :, 1]*(p[1:-1, 2:] - p_curr)*sn[:, :, 1, 0]/den[:, :, 1]
        uf[:, :, 2] = uf[:, :, 2] - df[:, :, 2]*(p[:-2, 1:-1] - p_curr)*sn[:, :, 2, 0]/den[:, :, 2]
        uf[:, :, 3] = uf[:, :, 3] - df[:, :, 3]*(p[1:-1, :-2] - p_curr)*sn[:, :, 3, 0]/den[:, :, 3]

        vf[:, :, 0] = vf[:, :, 0] - df[:, :, 0]*(p[2:, 1:-1] - p_curr)*sn[:, :, 0, 1]/den[:, :, 0]
        vf[:, :, 1] = vf[:, :, 1] - df[:, :, 1]*(p[1:-1, 2:] - p_curr)*sn[:, :, 1, 1]/den[:, :, 1]
        vf[:, :, 2] = vf[:, :, 2] - df[:, :, 2]*(p[:-2, 1:-1] - p_curr)*sn[:, :, 2, 1]/den[:, :, 2]
        vf[:, :, 3] = vf[:, :, 3] - df[:, :, 3]*(p[1:-1, :-2] - p_curr)*sn[:, :, 3, 1]/den[:, :, 3]

        return uf, vf, df

    def _correct(self, df):
        u = self.u[1:-1, 1:-1]
        v = self.v[1:-1, 1:-1]
        p = self.p
        p_corr = self.p_corr
        d = self.d[1:-1, 1:-1]
        uf, vf = self.uf, self.vf

        # Correct cell velocities
        p_corr_f = self._interpolate_faces(self.p_corr)
        dp_corr_x, dp_corr_y = self._compute_gradient(p_corr_f)

        if False:
            u[:, :] -= dp_corr_x*d
            v[:, :] -= dp_corr_y*d

        # Correct pressure
        if True:
            p[:, :] += self.p_corr*self.omega_p_corr

        # Finally must correct the face values
        sn, rc = self.sn, self.rc
        den = sn[:, :, :, 0]*rc[:, :, :, 0] + sn[:, :, :, 1]*rc[:, :, :, 1]
        p_curr = self.p_corr[1:-1, 1:-1]

        if False:

            uf[:, :, 0] -= uf[:, :, 0] - df[:, :, 0]*(p_corr[2:, 1:-1] - p_curr)*sn[:, :, 0, 0]/den[:, :, 0]
            uf[:, :, 1] -= uf[:, :, 1] - df[:, :, 1]*(p_corr[1:-1, 2:] - p_curr)*sn[:, :, 1, 0]/den[:, :, 1]
            uf[:, :, 2] -= uf[:, :, 2] - df[:, :, 2]*(p_curr - p_corr[:-2, 1:-1])*sn[:, :, 2, 0]/-den[:, :, 2]
            uf[:, :, 3] -= uf[:, :, 3] - df[:, :, 3]*(p_curr - p_corr[1:-1, :-2])*sn[:, :, 3, 0]/-den[:, :, 3]

            vf[:, :, 0] -= vf[:, :, 0] - df[:, :, 0]*(p_corr[2:, 1:-1] - p_curr)*sn[:, :, 0, 1]/den[:, :, 0]
            vf[:, :, 1] -= vf[:, :, 1] - df[:, :, 1]*(p_corr[1:-1, 2:] - p_curr)*sn[:, :, 1, 1]/den[:, :, 1]
            vf[:, :, 2] -= vf[:, :, 2] - df[:, :, 2]*(p_curr - p_corr[:-2, 1:-1])*sn[:, :, 2, 1]/-den[:, :, 2]
            vf[:, :, 3] -= vf[:, :, 3] - df[:, :, 3]*(p_curr - p_corr[1:-1, :-2])*sn[:, :, 3, 1]/-den[:, :, 3]
