import numpy as np

class Solver(object):
    def __init__(self, grid):
        self.grid = grid
        self.sn, self.rc = grid.core_face_norms, grid.core_links
        self.areas = grid.core_cell_areas

    def solve(self):
        raise NotImplementedError

    def _setup_bcs(self, bcs):
        raise NotImplementedError

    def _compute_interpolation_coeffs(self, **kwargs):
        grid = self.grid
        self.alpha, = grid.add_link_fields('alpha')
        alphas = self.alpha

        if kwargs.get('interpolation_method', 'volume_weighted') == 'volume_weighted':
            areas = grid.core_cell_areas

            alphas[:-1, :, 0] = areas[1:, :]/(areas[1:, :] + areas[:-1, :])
            alphas[:, :-1, 1] = areas[:, 1:]/(areas[:, 1:] + areas[:, :-1])
            alphas[1:, :, 2] = areas[:-1, :]/(areas[:-1, :] + areas[1:, :])
            alphas[:, 1:, 3] = areas[:, :-1]/(areas[:, :-1] + areas[:, 1:])
        else:
            raise ValueError

    def _interpolate_faces(self, cell_field):
        alphas = self.alpha
        face_field = np.ndarray(self.grid.link_shape, dtype=float, order='F')
        face_field[:, :, 0] = alphas[:, :, 0]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 0])*cell_field[2:, 1:-1]
        face_field[:, :, 1] = alphas[:, :, 1]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 1])*cell_field[1:-1, 2:]
        face_field[:, :, 2] = alphas[:, :, 2]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 2])*cell_field[:-2, 1:-1]
        face_field[:, :, 3] = alphas[:, :, 3]*cell_field[1:-1, 1:-1] + (1 - alphas[:, :, 3])*cell_field[1:-1, :-2]

        return face_field

    def _compute_gradient(self, face_field):
        sn = self.sn
        areas = self.areas

        grad_field_x = np.sum(face_field*sn[:, :, :, 0], dtype=float, axis=-1)/areas
        grad_field_y = np.sum(face_field*sn[:, :, :, 1], dtype=float, axis=-1)/areas

        return grad_field_x, grad_field_y
