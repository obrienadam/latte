from solvers.piso import Piso

class MultiphasePiso(Piso):
    def __init__(self, grid, **kwargs):
        super(MultiphasePiso, self).__init__(grid, **kwargs)

        self.f, self.k = grid.add_cell_fields('f', 'k')
        self.df_x, self.df_y = grid.add_cell_fields('df_x', 'df_y')

        self.rho_1 = kwargs['rho_1']
        self.rho_2 = kwargs['rho_2']
        self.mu_1 = kwargs['mu_1']
        self.mu_2 = kwargs['mu_2']

    def _advect_interface(self):
        pass

    def _compute_curvature(self):
        ff = self._interpolate_faces(self.f)
        self.df_x[1:-1, 1:-1], self.df_y[1:-1, 1:-1] = self._compute_gradient(ff)

    def _compute_rho(self):
        self.rho = self.f*self.rho_1 + (1 - self.f)*self.rho_2

    def _compute_mu(self):
        self.mu = self.rho/(self.f*self.rho_1/self.mu_1 + (1 - self.f)*self.rho_2/self.mu_2)