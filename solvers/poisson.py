from solver import Solver
from numpy import dot, zeros, transpose, array
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pylab as plt

class Poisson(Solver):
    def __init__(self, grid, **kwargs):
        super(Poisson, self).__init__(grid, **kwargs)

        grid.add_fields('phi', 'gamma', face_centered=True, cell_centered=True)

        # Constant fields
        grid.add_fields('d_e', 'd_w', 'd_n', 'd_s', 'd_p', 'b', cell_centered=True)

        # Get the cell face norms (these methods are super expensive!)
        sfe, sfw, sfn, sfs = grid.face_norms()

        # Get all of the relative cell vecs for interior cells
        rce, rcw, rcn, rcs = grid.cell_rvecs()

        nI, nJ = grid.interior_shape()
        a_e, a_w, a_n, a_s, a_p = grid.cell_data('d_e', 'd_w', 'd_n', 'd_s', 'd_p')

        for j in xrange(nJ):
            for i in xrange(nI):
                a_e[i,j] = dot(sfe[i,j], sfe[i,j])/dot(sfe[i,j], rce[i,j])
                a_w[i,j] = dot(sfw[i,j], sfw[i,j])/dot(sfw[i,j], rcw[i,j])
                a_n[i,j] = dot(sfn[i,j], sfn[i,j])/dot(sfn[i,j], rcn[i,j])
                a_s[i,j] = dot(sfs[i,j], sfs[i,j])/dot(sfs[i,j], rcs[i,j])

        for j in xrange(nJ):
            for i in xrange(nI):
                a_p[i,j] = -sum((a_e[i,j], a_w[i,j], a_n[i,j], a_s[i,j]))

        grid.set_const_field(gamma=kwargs.get('gamma', 1.))
        bcs = kwargs.get('bcs', None)

        if bcs:
            self.setup_bcs(bcs)

        print "Initialization of poisson solver complete."

    def setup_bcs(self, bcs):
        """
        :param bcs: Dictionary containing bc information
        :return:
        """
        a_e, a_w, a_n, a_s, a_p, b = self.grid.cell_data('d_e', 'd_w', 'd_n', 'd_s', 'd_p', 'b')

        if bcs['east']['type'] == 'zero_gradient':
            a_p[-1,:] += a_e[-1,:]
        elif bcs['east']['type'] == 'fixed':
            b[-1,:] -= a_e[-1,:]*bcs['east']['value']

        if bcs['west']['type'] == 'zero_gradient':
            a_p[0,:] += a_w[0,:]
        elif bcs['west']['type'] == 'fixed':
            b[0,:] -= a_w[0,:]*bcs['west']['value']

        if bcs['north']['type'] == 'zero_gradient':
            a_p[:,-1] += a_n[:,-1]
        elif bcs['north']['type'] == 'fixed':
            b[:,-1] -= a_n[:,-1]*bcs['north']['value']

        if bcs['south']['type'] == 'zero_gradient':
            a_p[:,0] += a_s[:,0]
        elif bcs['south']['type'] == 'fixed':
            b[:,0] -= a_s[:,0]*bcs['south']['value']

        a_e[-1,:] = 0.
        a_w[0,:] = 0.
        a_n[:,-1] = 0.
        a_s[:,0] = 0.


    def solve(self, **kwargs):
        a_e, a_w, a_n, a_s, a_p, b = self.grid.cell_data('d_e', 'd_w', 'd_n', 'd_s', 'd_p', 'b')
        phi = self.grid.cell_data('phi')
        nx, ny = self.grid.interior_shape()
        n = nx*ny

        d1 = a_s.reshape(n, order='F')[nx:]
        d2 = a_w.reshape(n, order='F')[1:]
        d3 = a_p.reshape(n, order='F')
        d4 = a_e.reshape(n, order='F')[0:n - 1]
        d5 = a_n.reshape(n, order='F')[0:n - nx]

        d2[nx - 1::nx] = 0.
        d4[nx - 1::nx] = 0.

        spmat = sp.diags([d1, d2, d3, d4, d5], [-nx, -1, 0, 1, nx], format='csr')

        rhs = b.reshape(n, order='F')

        phi = spla.spsolve(spmat, rhs).reshape((nx, ny), order='F')

        x, y = self.grid.cell_centers()

        plt.contourf(x, y, phi)
        plt.show()

        if kwargs.get('progress_bar', None):
            kwargs.get('progress_bar').setValue(100)




if __name__ == '__main__':
    from grid.finite_volume import FvGrid2D

    g = FvGrid2D((301, 301), (1, 1))


    bce = {'type': 'zero_gradient', 'value': 1.}
    bcw = {'type': 'fixed', 'value': 2.}
    bcn = {'type': 'zero_gradient', 'value': 3.}
    bcs = {'type': 'fixed', 'value': 1.}

    solver = Poisson(g, bcs={'east': bce, 'west':bcw, 'north':bcn, 'south': bcs}, gamma=0.01)
    solver.solve(progress_bar=None)