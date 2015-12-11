from poisson import Poisson
from numpy import dot
from grid.finite_volume import FvGrid2D

class Simple(Poisson):
    def __init__(self, grid, **kwargs):
        super(Simple, self).__init__(grid, **kwargs)

        grid.rename_field('phi', 'p')
        grid.add_fields('u', 'v', 'm', 'pCorr', 'rho', 'mu', cell_centered=True, face_centered=True)
        grid.add_fields('gradP_x', 'gradP_y', 'gradPCorr_x', 'gradPCorr_y', cell_centered=True)

        rho, mu = self.grid.cell_data('rho', 'mu')
        rho[:, :] = kwargs.get('rho', 1.)
        mu[:, :] = kwargs.get('mu', 0.1)

    def solve(self, **kwargs):
        progress_bar = kwargs.get('progress_bar', None)

        for iterNo in xrange(self.max_iters):

            compute_momentum(self.grid)
            correct_continuity(self.grid)

            if progress_bar:
                progress_bar.setValue(100.*(iterNo + 1)/self.max_iters)

def compute_momentum(grid):
    """
    Computes the momentum equation using a given pressure field
    :param grid: The computational grid
    :return:
    """
    # Get the mesh face vectors
    sfe, sfw, sfn, sfs = grid.east_face_norms(), grid.west_face_norms(), grid.north_face_norms(), grid.south_face_norms()
    # Get the relevant fields
    u, v, p, rho, mu = grid.cell_data('u', 'v', 'p', 'rho', 'mu')
    d_e, d_w, d_n, d_s, d_p, b = grid.cell_data('d_e', 'd_w', 'd_n', 'd_s', 'd_p', 'b')

    # Get the mass flow at the faces
    mh, mv = grid.face_data('m')

    num_i, num_j = grid.interior_shape()

    for j in xrange(num_j):
        for i in xrange(num_i):
            a_e, a_w, a_n, a_s = d_e[i, j], d_w[i, j], d_n[i, j], d_s[i, j]

            a_e +=  min(mv[i + 1, j], 0.)
            a_w += -max(mv[i, j], 0.)
            a_n +=  min(mh[i, j + 1], 0.)
            a_s += -max(mh[i, j], 0.)






def correct_continuity(grid):
    """
    Computes the pressure corrections and corrects the velocity field, mass flow field and pressure field
    :param grid: The computational grid
    :return:
    """
    u, v, m, p, pCorr, rho = grid.cell_data('u', 'v', 'm', 'p', 'pCorr', 'rho')

if __name__ == '__main__':
    grid = FvGrid2D((101, 101), (1, 1), num_ghost=2)
    solver = Simple(grid)
    solver.solve()