import numpy as np
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt

class FdGrid2D(object):
    def __init__(self, shape, dimensions):
        """
        Initialize a 2D finite difference mesh
        :param shape: The nodal resolution of the mesh
        :param dimensions: The absolution dimensions of the mesh
        :return: A constructed finite difference mesh
        """
        xdims = np.linspace(0, dimensions[0], shape[0])
        ydims = np.linspace(0, dimensions[1], shape[1])
        self.xnodes, self.ynodes = np.meshgrid(xdims, ydims, indexing='ij')

        self.shape = shape
        self.dimensions = dimensions
        self.fields = {}
        self.boundary_fields = {}

    def add_fields(self, *args):
        """
        Add fields to the finite difference mesh
        :param args: The names of the fields to be added
        :return:
        """
        for field_name in args:
            self.fields[field_name] = np.zeros(self.shape)
            field = self.fields[field_name]
            self.boundary_fields[field_name] = {'east': field[-1, :],
                                                'west': field[0, :],
                                                'north': field[:, -1],
                                                'south': field[:, 0]}

    def node_data(self, *args):
        """
        :param args: Names of fields
        :return: A tuple containing the node data of the desired fields
        """
        return tuple([self.fields[field_name] for field_name in args])

    def get_field_names(self):
        """
        :return: A list of the field names
        """
        return self.fields.keys()

    def node(self, i, j):
        """
        :param i: Index i
        :param j: Index j
        :return: A tuple containing the x, y coordinate of the node
        """
        return self.xnodes[i, j], self.ynodes[i, j]

    def nodes(self):
        """
        :return: A tuple containing all of the node coordinates
        """
        return self.xnodes, self.ynodes

    def plot(self, **kwargs):
        plt.plot(self.xnodes, self.ynodes,
                 np.transpose(self.xnodes), np.transpose(self.ynodes),
                 color='gray', linewidth=1.5)

        if kwargs.get('mark_nodes', False):
            plt.plot(self.xnodes, self.ynodes,
                     's', color='black', markersize=3)

        minx = np.min(np.min(self.xnodes))
        maxx = np.max(np.max(self.xnodes))
        miny = np.min(np.min(self.ynodes))
        maxy = np.max(np.max(self.ynodes))

        minx -= 0.1*self.dimensions[0]
        maxx += 0.1*self.dimensions[0]
        miny -= 0.1*self.dimensions[1]
        maxy += 0.1*self.dimensions[1]

        plt.axis((minx, maxx, miny, maxy))
        plt.axes().set_aspect('equal')

    def plot_show(self):
        plt.show()

if __name__ == '__main__':
    grid = FdGrid2D((21, 21), (2, 1))
    grid.add_fields('u', 'v', 'p', 'm', 'pCorr')

    grid.plot(mark_nodes=True)
    grid.plot_show()