from solvers.simple import Simple

class Piso(Simple):
    def __init__(self, grid, **kwargs):
        super(Piso, self).__init__(grid, **kwargs)
        self._num_p_corrs = kwargs.get('num_p_corrs', 2)

    def solve(self):
        self._solve_momentum()

        for i in xrange(self._num_p_corrs):
            self._solve_p_corr()