import numpy as np
from geometry import funcs as gf

class Field(np.ndarray):
    def __init__(self, shape, val=0.):
        super(Field, self).__init__(shape, dtype=float, order='F')
        self.fill(val)

    @property
    def core(self):
        return self[1:-1, 1:-1]

class FvQuadrilateralGrid(object):
    """
    This is the generalized finite volume quadrilateral grid, which can be initialized simply by passing it a 2D array
    of nodes representing cell corners.
    """
    array_config = {'dtype': float, 'order': 'F'}

    def __init__(self, corner_nodes):
        self._shape = np.array(corner_nodes[:, :, 0].shape) + 1

        self._corner_nodes = corner_nodes

        self._cell_nodes = np.zeros((self.shape[0], self.shape[1], 2), **self.array_config)
        self._cell_areas = np.zeros(self.shape, **self.array_config)
        self._face_norms = np.zeros((self.shape[0], self.shape[1], 4, 2), **self.array_config)

        # Initialize node properties
        self._initialize_core_cells(corner_nodes)
        self._initialize_boundary_cells(corner_nodes)

        # Initialize links
        self._initialize_links()

        self._cell_fields = {'cell_areas': self._cell_areas}
        self._link_fields = {}

    def _initialize_core_cells(self, corner_nodes):
        cell_nodes = self._cell_nodes[1:-1, 1:-1]
        cell_areas = self._cell_areas[1:-1, 1:-1]
        face_norms = self._face_norms[1:-1, 1:-1]

        num_i, num_j = corner_nodes[:, :, 0].shape
        nodes_x, nodes_y = corner_nodes[:, :, 0], corner_nodes[:, :, 1]
        for j in xrange(num_j - 1):
            for i in xrange(num_i - 1):
                cell_polygon = (nodes_x[i, j], nodes_x[i + 1, j], nodes_x[i + 1, j + 1], nodes_x[i, j + 1]), \
                               (nodes_y[i, j], nodes_y[i + 1, j], nodes_y[i + 1, j + 1], nodes_y[i, j + 1])

                cell_nodes[i, j] = np.array(gf.polygon_centroid(*cell_polygon))
                cell_areas[i, j] = gf.polygon_area(*cell_polygon)
                face_norms[i, j] = np.array(gf.polygon_face_normals(*cell_polygon, start_point=1))

    def _initialize_boundary_cells(self, corner_nodes):
        nodes_x, nodes_y = corner_nodes[:, :, 0], corner_nodes[:, :, 1]
        self._cell_nodes[0, 0] = np.array((nodes_x[0, 0], nodes_y[0, 0]))
        self._cell_nodes[-1, 0] = np.array((nodes_x[-1, 0], nodes_y[-1, 0]))
        self._cell_nodes[-1, -1] = np.array((nodes_x[-1, -1], nodes_y[-1, -1]))
        self._cell_nodes[0, -1] = np.array((nodes_x[0, -1], nodes_y[0, -1]))

        right_cell_nodes = self._cell_nodes[-1, 1:-1]
        left_cell_nodes = self._cell_nodes[0, 1:-1]
        top_cell_nodes = self._cell_nodes[1:-1, -1]
        bottom_cell_nodes = self._cell_nodes[1:-1, 0]

        right_cell_nodes[:, 0] = nodes_x[-1, :-1] + 0.5*np.diff(nodes_x[-1, :])
        right_cell_nodes[:, 1] = nodes_y[-1, :-1] + 0.5*np.diff(nodes_y[-1, :])

        left_cell_nodes[:, 0] = nodes_x[0, :-1] + 0.5*np.diff(nodes_x[0, :])
        left_cell_nodes[:, 1] = nodes_y[0, :-1] + 0.5*np.diff(nodes_y[0, :])

        top_cell_nodes[:, 0] = nodes_x[:-1, -1] + 0.5*np.diff(nodes_x[:, -1])
        top_cell_nodes[:, 1] = nodes_y[:-1, -1] + 0.5*np.diff(nodes_y[:, -1])

        bottom_cell_nodes[:, 0] = nodes_x[:-1, 0   ] + 0.5*np.diff(nodes_x[:, 0])
        bottom_cell_nodes[:, 1] = nodes_y[:-1, 0] + 0.5*np.diff(nodes_y[:, 0])

    def _initialize_links(self):
        cell_nodes = self._cell_nodes

        self._link_vecs = np.zeros((self.shape[0] - 2, self.shape[1] - 2, 4, 2), **self.array_config)
        for j in xrange(1, cell_nodes.shape[1] - 1):
            for i in xrange(1, cell_nodes.shape[0] - 1):
                ic, jc = i - 1, j - 1

                self._link_vecs[ic, jc, 0] = cell_nodes[i + 1, j] - cell_nodes[i, j]
                self._link_vecs[ic, jc, 1] = cell_nodes[i, j + 1] - cell_nodes[i, j]
                self._link_vecs[ic, jc, 2] = cell_nodes[i - 1, j] - cell_nodes[i, j]
                self._link_vecs[ic, jc, 3] = cell_nodes[i, j - 1] - cell_nodes[i, j]


    def add_cell_fields(self, *args):
        new_fields = []

        for field_name in args:
            if field_name not in self._cell_fields:
                self._cell_fields[field_name] = np.zeros(self.shape, **self.array_config)

            new_fields.append(self._cell_fields[field_name])

        return tuple(new_fields)

    def add_link_fields(self, *args):
        new_fields = []

        for field_name in args:
            if field_name not in self._link_fields:
                self._link_fields[field_name] = np.zeros(self.link_shape, **self.array_config)

            new_fields.append(self._link_fields[field_name])

        return tuple(new_fields)

    def add_existing_cell_fields(self, *args):
        for field_tuple in args:
            field_name = field_tuple[0]
            field = field_tuple[1]

            if np.array_equal(field.shape, self.shape):
                self._cell_fields[field_name] = field
            elif np.array_equal(field.shape, self.core_shape):
                new_field, = self.add_cell_fields(field_name)
                new_field[1:-1, 1:-1] = field
            else:
                raise AttributeError

    def get_cell_fields(self, *args):
        fields = []

        for field_name in args:
            fields.append(self._cell_fields[field_name])

        return tuple(fields)

    def get_core_cell_fields(self, *args):
        fields = []

        for field_name in args:
            fields.append(self._fields[field_name][1:-1, 1:-1])

        return tuple(fields)

    def get_link_fields(self, *args):
        fields = []

        for field_name in args:
            fields.append(self._link_fields[field_name])

        return tuple(fields)

    @property
    def shape(self):
        return self._shape

    @property
    def core_shape(self):
        return self._shape[0] - 2, self._shape[1] - 2

    @property
    def link_shape(self):
        return self.core_shape + (4,)

    @property
    def cell_nodes(self):
        return self._cell_nodes

    @property
    def corner_nodes(self):
        return self._corner_nodes

    @property
    def core_cell_nodes(self):
        return self._cell_nodes[1:-1, 1:-1]

    @property
    def core_cell_areas(self):
        return self._cell_areas[1:-1, 1:-1]

    @property
    def core_face_norms(self):
        return self._face_norms[1:-1, 1:-1]

    @property
    def core_links(self):
        return self._link_vecs

class FvRectilinearGrid(FvQuadrilateralGrid):
    def __init__(self, shape, dims):
        x, y = np.linspace(0., dims[0], shape[0] + 1), np.linspace(0., dims[1], shape[1] + 1)
        x, y = np.meshgrid(x, y, indexing='ij')
        corner_nodes = np.dstack((x, y))

        super(FvRectilinearGrid, self).__init__(corner_nodes)

class FvEquidistantGrid(FvQuadrilateralGrid):
    """
    This class is a 2D cartesian grid with equal spacing. It has the advantage that it can initialize much faster than
    the generalized quadrilateral grid.
    """
    def __init__(self, n, width):
        self._shape = np.array((n, n)) + 2

        self.width = width
        self.h = float(width)/n

        self._cell_nodes = np.zeros((self.shape[0], self.shape[1], 2), **self.array_config)
        self._cell_areas = np.zeros(self.shape, **self.array_config)
        self._face_norms = np.zeros((self.shape[0], self.shape[1], 4, 2), **self.array_config)

        nodes = np.linspace(0, width, n + 1)
        nodes_x, nodes_y = np.meshgrid(nodes, nodes, indexing='ij')

        self._initialize_core_cells(nodes_x, nodes_y)
        self._initialize_boundary_cells(nodes_x, nodes_y)
        self._initialize_links()

        self._cell_fields = {'cell_areas': self._cell_areas}
        self._link_fields = {}

    def _initialize_core_cells(self, nodes_x, nodes_y):
        cell_nodes = self._cell_nodes[1:-1, 1:-1]
        cell_areas = self._cell_areas[1:-1, 1:-1]
        face_norms = self._face_norms[1:-1, 1:-1]

        cell_nodes[:, :, 0] = nodes_x[:-1, :-1] + 0.5*self.h
        cell_nodes[:, :, 1] = nodes_y[:-1, :-1] + 0.5*self.h

        cell_areas[:, :] = self.h**2

        face_norms[:, :, 0] = np.array((self.h, 0.))
        face_norms[:, :, 1] = np.array((0., self.h))
        face_norms[:, :, 2] = np.array((-self.h, 0.))
        face_norms[:, :, 3] = np.array((0., -self.h))

    def _initialize_links(self):
        cell_nodes = self._cell_nodes

        self._link_vecs = np.zeros((self.shape[0] - 2, self.shape[1] - 2, 4, 2), **self.array_config)

        self._link_vecs[:, :, 0] = cell_nodes[2:, 1:-1] - cell_nodes[1:-1, 1:-1]
        self._link_vecs[:, :, 1] = cell_nodes[1:-1, 2:] - cell_nodes[1:-1, 1:-1]
        self._link_vecs[:, :, 2] = cell_nodes[:-2, 1:-1] - cell_nodes[1:-1, 1:-1]
        self._link_vecs[:, :, 3] = cell_nodes[1:-1, :-2] - cell_nodes[1:-1, 1:-1]

if __name__ == '__main__':
    field = Field((10, 10))
    print field.core