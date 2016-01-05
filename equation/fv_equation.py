import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

class Term(object):
    def __init__(self, shape):
        self.__shape = shape
        self.a = np.zeros(self.shape + (5,), dtype=float, order='F')
        self.b = np.zeros(self.shape, dtype=float, order='F')

    @property
    def shape(self):
        return self.__shape

    def __add__(self, other):
        term = Term(self.shape)
        term.a = self.a + other.a
        term.b = self.b + other.b

        return term

    def __sub__(self, other):
        term = Term(self.shape)
        term.a = self.a - other.a
        term.b = self.b - other.b

        return term

class DiffusionTerm(Term):
    def __init__(self, grid, *args):
        super(DiffusionTerm, self).__init__(grid.core_shape)

        sn, rc = grid.core_face_norms, grid.core_links

        self.gamma = np.ones(self.shape + (4,), dtype=float, order='F')

        for coeff in args:
            self.gamma *= coeff

        for i in xrange(4):
            num = sn[:, :, i, 0]*sn[:, :, i, 0] + sn[:, :, i, 1]*sn[:, :, i, 1]
            den = sn[:, :, i, 0]*rc[:, :, i, 0] + sn[:, :, i, 1]*rc[:, :, i, 1]
            self.a[:, :, i] = self.gamma[:, :, i]*num/den

        self.a[:, :, 4] = -np.sum(self.a[:, :, :-1], axis=2, dtype=float)

class AdvectionTerm(Term):
    def __init__(self, grid, u, v, *args, **kwargs):
        super(AdvectionTerm, self).__init__(grid.core_shape)

        sn, rc = grid.core_face_norms, grid.core_links

        isnum = lambda num: isinstance(num, float) or isinstance(num, int)

        if isnum(u):
            u = u*np.ones(self.shape + (4,), dtype=float, order='F')
        if isnum(v):
            v = v*np.ones(self.shape + (4,), dtype=float, order='F')

        self.u = u
        self.v = v

        for coeff in args:
            self.u *= coeff
            self.v *= coeff

        scheme = kwargs.get('scheme', 'central')

        if scheme is 'central':
            for i in xrange(4):
                self.a[:, :, i] = 0.5*(sn[:, :, i, 0]*self.u[:, :, i] + sn[:, :, i, 1]*self.v[:, :, i])

            self.a[:, :, 4] = np.sum(self.a[:, :, :-1], axis=2, dtype=float)

        elif scheme is 'upwind':
            upwind_nb = np.vectorize(lambda a: min(a, 0.))
            upwind_p = np.vectorize(lambda a: max(a, 0.))

            self.a[:, :, 4] = 0.

            for i in xrange(4):
                dot_us = sn[:, :, i, 0]*self.u[:, :, i] + sn[:, :, i, 1]*self.v[:, :, i]
                self.a[:, :, i] = upwind_nb(dot_us)
                self.a[:, :, 4] += upwind_p(dot_us)

        else:
            raise ValueError

class TemporalTerm(Term):
    def __init__(self, shape):
        super(TemporalTerm, self).__init__(shape)

class FvEquation(object):
    def __init__(self, var, **kwargs):
        self.var = var
        self.shape = var.shape

        self.bcs = kwargs.get('bcs', {'type': ['fixed']*4,
                                      'value': np.array(range(0, 4)),
                                      })

        # Implicit coefficients
        self.a = np.zeros((self.shape[0], self.shape[1], 5), dtype=float)
        self.a_core = self.a[1:-1, 1:-1]
        self.a_boundary = self.a[-1, :, :], self.a[:, -1, :], self.a[0, :, :], self.a[:, 0, :]
        self.a_nb = [2, 3, 0, 1]

        # Source term
        self.b = np.zeros(self.shape, dtype=float)
        self.b_core = self.b[1:-1, 1:-1]
        self.b_boundary = self.b[-1, :], self.b[:, -1], self.b[0, :], self.b[:, 0]

    def solve(self):
        # Apply bcs
        bc_types = self.bcs['type']
        bc_vals = self.bcs['value']

        for bc_coeff, bc_source, bc_type, bc_value, nb_coeff  in zip(self.a_boundary, self.b_boundary, bc_types, bc_vals, self.a_nb):
            if bc_type is 'fixed' or 'f':
                bc_coeff[:, 4] = 1.
            elif bc_type is 'normal_gradient' or 'ng':
                bc_coeff[:, 4] = 1.
                bc_coeff[:, nb_coeff] = -1.
            else:
                raise ValueError

            bc_source[:] = bc_value

        num_i, num_j = self.var.shape
        num = num_i*num_j

        diags = self.a[:, :, 0].flatten(order='F')[:-1], \
                self.a[:, :, 1].flatten(order='F')[0:num - num_i], \
                self.a[:, :, 2].flatten(order='F')[1:], \
                self.a[:, :, 3].flatten(order='F')[num_i:], \
                self.a[:, :, 4].flatten(order='F')

        spmat = sp.diags(diags, [1, num_i, -1, -num_i, 0], format='csr')
        rhs = self.b.flatten(order='F')

        self.var[:, :] = spla.spsolve(spmat, rhs).reshape(self.var.shape, order='F')

        return self.var

    def __eq__(self, term):
        self.a_core[:, :] = term.a
        self.b_core[:, :] = term.b