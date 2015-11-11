import numpy as np
from finite_difference import FdGrid2D
import matplotlib.pyplot as plt

class FvGrid2D(FdGrid2D):
    def __init__(self, shape, dimensions):
        """
        Initialize a 2D finite volume grid
        :param shape: The nodal resolution of the mesh
        :param dimensions: The absolute dimensions of the mesh
        :return: A constructed finite volume mesh
        """
        super(FvGrid2D, self).__init__(shape, dimensions)

        xdims = np.linspace(0, dimensions[0], shape[0])
        ydims = np.linspace(0, dimensions[1], shape[1])
        deltaxs = np.diff(xdims)
        deltays = np.diff(ydims)

        xcdims = [xdim + 0.5*deltax for xdim, deltax in zip(xdims, deltaxs)]
        ycdims = [ydim + 0.5*deltay for ydim, deltay in zip(ydims, deltays)]

        # Compute cell and face centroids
        self.cell_centers_x, self.cell_centers_y = np.meshgrid(xcdims, ycdims, indexing='ij')
        self.hface_centers_x, self.hface_centers_y = np.meshgrid(xcdims, ydims, indexing='ij')
        self.vface_centers_x, self.vface_centers_y = np.meshgrid(xdims, ycdims, indexing='ij')

        # Compute cell and face relative displacements
        self.rcell_h_x = np.diff(self.cell_centers_x, axis=0)
        self.rcell_h_y = np.diff(self.cell_centers_y, axis=0)
        self.rcell_v_x = np.diff(self.cell_centers_x, axis=1)
        self.rcell_v_y = np.diff(self.cell_centers_y, axis=1)

        self.rface_east_x = self.vface_centers_x[1:, :] - self.cell_centers_x
        self.rface_east_y = self.vface_centers_y[1:, :] - self.cell_centers_y

        self.rface_west_x = self.vface_centers_x[0:-1, :] - self.cell_centers_x
        self.rface_west_y = self.vface_centers_y[0:-1, :] - self.cell_centers_y

        self.rface_north_x = self.hface_centers_x[:, 1:] - self.cell_centers_x
        self.rface_north_y = self.hface_centers_y[:, 1:] - self.cell_centers_y

        self.rface_south_x = self.hface_centers_x[:, 0:-1] - self.cell_centers_x
        self.rface_south_y = self.hface_centers_y[:, 0:-1] - self.cell_centers_y

        # Compute the face parameters
        self.hface_norms_x = np.diff(self.ynodes, axis=0)
        self.hface_norms_y = -np.diff(self.xnodes, axis=0)

        self.vface_norms_x = np.diff(self.ynodes, axis=1)
        self.vface_norms_y = -np.diff(self.xnodes, axis=1)

        # Compute diffusion metrics

        self.cell_shape = len(xcdims), len(ycdims)
        self.hface_shape = len(xcdims), len(ydims)
        self.vface_shape = len(xdims), len(ycdims)

    def add_fields(self, *args, **kwargs):
        """
        Add fields to the finite volume mesh
        :param args: The names of the fields to be added
        :param kwargs: Locations of the fields to be added
        :return:
        """
        for field_name in args:
            field = self.fields[field_name] = {}

            if kwargs.get('face_centered', False):
                kwargs['hface_centered'] = True
                kwargs['vface_centered'] = True

            field['cell_data'] = np.zeros(self.cell_shape) if kwargs.get('cell_centered', False) else None
            field['node_data'] = np.zeros(self.shape) if kwargs.get('node_centered', False) else None
            field['hface_data'] = np.zeros(self.hface_shape) if kwargs.get('hface_centered', False) else None
            field['vface_data'] = np.zeros(self.vface_shape) if kwargs.get('vface_centered', False) else None

            if field['cell_data'] != None:
                cell_data = field['cell_data']
                field['boundary_data'] = {'east': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[-1, :]},
                                            'west': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[0, :]},
                                            'north': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[:, -1]},
                                            'south': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[:, 0]}}

    def cell_data(self, *args):
        """
        :param args: Names of fields
        :return: A tuple containing the cell data of the desired fields
        """
        return tuple([self.fields[field_name]['cell_data'] for field_name in args])

    def node_data(self, *args):
        """
        :param args: Names of fields
        :return: A tuple containing the node data of the desired fields
        """
        return tuple([self.fields[field_name]['node_data'] for field_name in args])

    def face_data(self, field_name):
        """
        :param args: Names of fields
        :return: A tuple containing tuples with the horizontal and vertical face data of the desired fields
        """
        return self.fields[field_name]['hface_data'], self.fields[field_name]['vface_data']

    def xc(self, i, j):
        """
        :param i: Index i
        :param j: Index j
        :return: Tuple containing the x, y coordinate of the cell center
        """
        return self.cell_centers_x[i, j], self.cell_centers_y[i, j]

    def sfe(self, i, j):
        """
        Get the east face normal at cell i, j
        :param i: Index i
        :param j: Index j
        :return: Tuple containing x, y components of the outward normal
        """
        return self.vface_norms_x[i + 1, j], self.vface_norms_x[i + 1, j]

    def sfw(self, i, j):
        """
        Get the west face normal at cell i, j
        :param i: Index i
        :param j: Index j
        :return: Tuple containing x, y components of the outward normal
        """
        return -self.vface_norms_x[i, j], -self.vface_norms_x[i, j]

    def sfn(self, i, j):
        """
        Get the north face normal at cell i, j
        :param i: Index i
        :param j: Index j
        :return: Tuple containing x, y components of the outward normal
        """
        return self.hface_norms_x[i, j + 1], self.vface_norms_y[i, j + 1]

    def sfs(self, i, j):
        """
        Get the south face normal at cell i, j
        :param i: Index i
        :param j: Index j
        :return: Tuple containing x, y components of the outward normal
        """
        return -self.hface_norms_x[i, j], -self.hface_norms_x[i, j]

    def east_face_norms(self):
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.vface_norms_x[1:, :], self.vface_norms_y[1:, :])])

    def west_face_norms(self):
        return [zip(xvals, yvals) for xvals, yvals in zip(-self.vface_norms_x[:-1, :], -self.vface_norms_y[:-1, :])]

    def north_face_norms(self):
        return [zip(xvals, yvals) for xvals, yvals in zip(self.hface_norms_x[:, 1:], self.hface_norms_y[:, 1:])]

    def south_face_norms(self):
        return [zip(xvals, yvals) for xvals, yvals in zip(-self.hface_norms_x[:, :-1], -self.hface_norms_y[:, :-1])]

    def cell_centers(self):
        """
        :return: A tuple containing all the coordinates of the cell centers
        """
        return self.cell_centers_x, self.cell_centers_y

    def plot(self, **kwargs):
        super(FvGrid2D, self).plot(**kwargs)

        if kwargs.get('mark_cell_centers', False):
            plt.plot(self.cell_centers_x, self.cell_centers_y,
                     'o', color='red', markersize=3)

        if kwargs.get('mark_faces', False):
            plt.plot(self.hface_centers_x, self.hface_centers_y,
                     'p', color='blue', markersize=3)
            plt.plot(self.vface_centers_x, self.vface_centers_y,
                     'p', color='blue', markersize=3)

if __name__ == '__main__':
    grid = FvGrid2D((21, 21), (1, 1))
    grid.add_fields('u', 'v', 'p', 'pCorr', 'm', 'phi', cell_centered=True, face_centered=True)

    sfe = grid.east_face_norms()
    sfw = grid.west_face_norms()

    grid.plot(mark_nodes=True, mark_cell_centers=True, mark_faces=True)
    grid.plot_show()
