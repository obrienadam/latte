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

        self.x_cell_centers, self.y_cell_centers = np.meshgrid(xcdims, ycdims, indexing='ij')
        self.x_hface_centers, self.y_hface_centers = np.meshgrid(xcdims, ydims, indexing='ij')
        self.x_vface_centers, self.y_vface_centers = np.meshgrid(xdims, ycdims, indexing='ij')

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

    def face_data(self, args):
        """
        :param args: Names of fields
        :return: A tuple containing tuples with the horizontal and vertical face data of the desired fields
        """
        return tuple([(self.fields[field_name]['hface_data'], self.fields[field_name]['vface_data']) for field_name in args])

    def cell_center(self, index):
        """
        :param index: The (i, j) index of the cell
        :return: A tuple containing the coordinates of the cell center
        """
        return self.x_cell_centers[index], self.y_cell_centers[index]

    def cell_centers(self):
        """
        :return: A tuple containing all the coordinates of the cell centers
        """
        return self.x_cell_centers, self.y_cell_centers

    def plot(self, **kwargs):
        super(FvGrid2D, self).plot(**kwargs)

        if kwargs.get('mark_cell_centers', False):
            plt.plot(self.x_cell_centers, self.y_cell_centers,
                     'o', color='red', markersize=3)

        if kwargs.get('mark_faces', False):
            plt.plot(self.x_hface_centers, self.y_hface_centers,
                     'p', color='blue', markersize=3)
            plt.plot(self.x_vface_centers, self.y_vface_centers,
                     'p', color='blue', markersize=3)

if __name__ == '__main__':
    grid = FvGrid2D((21, 21), (1, 1))
    grid.add_fields('u', 'v', 'p', 'pCorr', 'm', 'phi', cell_centered=True, face_centered=True)
    grid.add_fields('')
    grid.plot(mark_nodes=True, mark_cell_centers=True, mark_faces=True)
    grid.plot_show()
