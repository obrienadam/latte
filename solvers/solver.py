class Solver(object):
    def __init__(self, grid, **kwargs):
        self.grid = grid

        self.max_iters = kwargs.get('max_iters', 1)
        self.time_accurate = kwargs.get('time_accurate', False)

    def solve(self, progress_bar):
        raise NotImplementedError
