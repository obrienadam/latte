from solver import Solver
from numpy import dot, array, reshape
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib.pylab as plt

class Poisson(Solver):
    def __init__(self, grid, **kwargs):
        super(Poisson, self).__init__(grid, **kwargs)

        grid.add_fields('phi', face_centered=True, cell_centered=True)

        # Constant fields
        grid.add_fields('a_e', 'a_w', 'a_n', 'a_s', 'a_p', cell_centered=True)
        a_e, a_w, a_n, a_s, a_p = grid.cell_data('a_e', 'a_w', 'a_n', 'a_s', 'a_p')

        # Get the cell face norms, however we don't want the last one!
        sfe = grid.east_face_norms()
        sfw = grid.west_face_norms()
        sfn = grid.north_face_norms()
        sfs = grid.south_face_norms()

        # Get all of the relative cell vecs for interior cells
        rce = grid.east_cell_rvecs()
        rcw = grid.west_cell_rvecs()
        rcn = grid.north_cell_rvecs()
        rcs = grid.south_cell_rvecs()
        nI, nJ = grid.cell_shape

        for j in xrange(nJ):
            for i in xrange(nI):
                a_e[i,j] = dot(sfe[i,j], sfe[i,j])/dot(sfe[i,j], rce[i,j]) if i < nI - 1 else 0.
                a_w[i,j] = dot(sfw[i,j], sfw[i,j])/dot(sfw[i,j], rcw[i-1,j]) if i > 0 else 0.
                a_n[i,j] = dot(sfn[i,j], sfn[i,j])/dot(sfn[i,j], rcn[i,j]) if j < nJ - 1 else 0.
                a_s[i,j] = dot(sfs[i,j], sfs[i,j])/dot(sfs[i,j], rcs[i,j-1]) if j > 0 else 0.

        for j in xrange(nJ):
            for i in xrange(nI):
                a_p[i,j] = -sum((a_e[i,j], a_w[i,j], a_n[i,j], a_s[i,j]))

        print "Initialization of poisson solver complete."


    def solve(self, progress_bar):
        a_e, a_w, a_n, a_s, a_p = self.grid.cell_data('a_e', 'a_w', 'a_n', 'a_s', 'a_p')
        nx, ny = self.grid.cell_shape
        n = nx*ny

        d1 = a_s.reshape(n, order='C')[nx:]
        d2 = a_w.reshape(n, order='C')[1:]
        d3 = a_p.reshape(n, order='C')
        d4 = a_e.reshape(n, order='C')[1:]
        d5 = a_n.reshape(n, order='C')[nx:]

        spmat = sp.diags([d1, d2, d3, d4, d5], [-nx, -1, 0, 1, nx], format='csr')
        plt.spy(spmat)
        plt.show()




if __name__ == '__main__':
    from grid.finite_volume import FvGrid2D

    g = FvGrid2D((11, 11), (1, 1))

    solver = Poisson(g)
    solver.solve(3)