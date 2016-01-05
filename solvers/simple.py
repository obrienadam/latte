from equation.fv_equation import FvEquation, DiffusionTerm, AdvectionTerm

class Simple(object):
    def __init__(self, grid, **kwargs):
        self.grid = grid

        self.u, self.v, self.p, self.p_corr = grid.add_fields('u', 'v', 'p', 'p_corr', cell_centered=True, link_centered=True)

        self.u_eqn = FvEquation(self.u)
        self.v_eqn = FvEquation(self.v)
        self.p_corr_eqn = FvEquation(self.p_corr)

    def _setup_bcs(self, bcs):
        self.bcs = bcs


