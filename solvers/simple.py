from poisson import Poisson
from numpy import dot
from grid.finite_volume import FvGrid2D

class Simple(Poisson):
    def __init__(self, grid, **kwargs):
        super(Simple, self).__init__(grid, **kwargs)

        self.omega_momentum = kwargs.get('omega_momentum', 1.)
        self.omega_pcorr = kwargs.get('omega_pcorr', 1.)
        self.bcs = kwargs.get('bcs')

        grid.rename_field(phi='p', gamma='mu')
        grid.add_fields('u', 'v', 'm', 'pCorr', 'rho', 'mu', cell_centered=True, face_centered=True)
        grid.add_fields('gradP_x', 'gradP_y', 'gradPCorr_x', 'gradPCorr_y', cell_centered=True)
        grid.add_fields('a_e', 'a_w', 'a_n', 'a_s', 'a_p', cell_centered=True)

        # Store relevant grid fields for easier and faster reference
        self.primitives = grid.cell_data('u', 'v', 'p')
        self.diff_coeffs = grid.cell_data('d_e', 'd_w', 'd_n', 'd_s', 'd_p')
        self.implicit_coeffs = grid.cell_data('a_e', 'a_w', 'a_n', 'a_s', 'a_p', 'b')
        self.mass_flux = grid.face_data('m')
        self.face_coeffs = {'rho': grid.face_data('rho'), 'mu': grid.face_data('mu')}

        # Store the most frequently used mesh metrics for convenience
        self.face_norms = grid.face_norms()
        self.cell_rvecs = grid.cell_rvecs()

        grid.set_const_field(rho=kwargs.get('rho', 1.), mu=kwargs.get('mu', 0.1))

    def solve(self, **kwargs):
        progress_bar = kwargs.get('progress_bar', None)

        for iterNo in xrange(self.max_iters):
            self.compute_momentum()

            if progress_bar:
                progress_bar.setValue(100.*(iterNo + 1)/self.max_iters)

    def compute_momentum(self):
        """
        Computes the momentum equation using a given pressure field
        :param grid: The computational grid
        :return:
        """
        grid = self.grid
        sfe, sfw, sfn, sfs = self.face_norms
        rce, rcw, rcn, rcs = self.cell_rvecs

        # Get the relevant fields
        u, v, p = self.primitives
        d_e, d_w, d_n, d_s, d_p = self.diff_coeffs
        a_e, a_w, a_n, a_s, a_p, b = self.implicit_coeffs

        mu_e, mu_w, mu_n, mu_s = self.face_coeffs['mu']

        # Get the mass flow at the faces
        m_e, m_w, m_n, m_s = self.mass_flux

        num_i, num_j = grid.interior_shape()
        for j in xrange(num_j):
            for i in xrange(num_i):
                a_e[i, j], a_w[i, j], a_n[i, j], a_s[i, j] = mu_e[i, j]*d_e[i, j], \
                                                             mu_w[i, j]*d_w[i, j], \
                                                             mu_n[i, j]*d_n[i, j], \
                                                             mu_s[i, j]*d_s[i, j]

                a_e[i, j] +=  min(m_e[i, j], 0.)
                a_w[i, j] += -max(m_w[i, j], 0.)
                a_n[i, j] +=  min(m_n[i, j], 0.)
                a_s[i, j] += -max(m_s[i, j], 0.)

                a_p[i, j] = max(m_e[i, j], 0.) + mu_e[i, j]*d_e[i, j] \
                            - min(m_w[i, j], 0.) + mu_w[i, j]*d_w[i, j] \
                            + max(m_n[i, j], 0.) + mu_n[i, j]*d_n[i, j] \
                            - min(m_s[i, j], 0.) + mu_s[i, j]*d_s[i, j]
                a_p[i, j] /= self.omega_momentum

        self.set_bcs(self.bcs)

        self.rhie_chow_interpolation()

    def setup_bcs(self, bcs):
        a_e, a_w, a_n, a_s, a_p, b = self.implicit_coeffs
        types = bcs['east']['type'],

        for coeff, type in zip((a_e[-1, :], a_w[0, :], a_n[:, -1], a_s[:, 0]), types):
            pass

    def rhie_chow_interpolation(self):
        pass

if __name__ == '__main__':
    grid = FvGrid2D((201, 201), (1, 1), num_ghost=2)
    solver = Simple(grid, omega_momentum=0.5, omega_pcorr=0.5)
    solver.solve()