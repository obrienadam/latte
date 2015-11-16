from .poisson import Poisson
from numpy import dot

class Simple(Poisson):
    def __init__(self, grid, **kwargs):
        super(Simple, self).__init__(grid, **kwargs)

        grid.rename_field('phi', 'p')
        grid.add_fields('u', 'v', 'm', 'pCorr', cell_centered=True, face_centered=True)
        grid.add_fields('gradP_x', 'gradP_y', 'gradPCorr_x', 'gradPCorr_y', cell_centered=True)

        self.rho = kwargs.get('rho', 1.)
        self.mu = kwargs.get('mu', 0.1)

    def solve(self, progress_bar):

        for iterNo in xrange(self.max_iters):

            compute_momentum(self.grid)
            correct_continuity(self.grid)

            progress_bar.setValue(100.*(iterNo + 1)/self.max_iters)

def compute_momentum(grid):
    """
    Computes the momentum equation using a given pressure field
    :param grid: The computational grid
    :return:
    """
    sfe = grid.east_face_norms()
    sfw = grid.west_face_norms()
    sfn = grid.north_face_norms()
    sfs = grid.south_face_norms()

    u, v, p = grid.cell_data('u', 'v', 'p')

    num_i, num_j = grid.cell_shape

    for j in xrange(num_j):
        for i in xrange(num_i):
            pass



def correct_continuity(grid):
    """
    Computes the pressure corrections and corrects the velocity field, mass flow field and pressure field
    :param grid: The computational grid
    :return:
    """
    u, v, m, p, pCorr = grid.cell_data('u', 'v', 'm', 'p', 'pCorr')


