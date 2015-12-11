import numpy as np
import matplotlib.pyplot as plt

class FdGrid2D(object):
    def __init__(self, shape, dimensions, num_ghost=1):
        """
        Initialize a 2D finite difference mesh
        :param shape: The nodal resolution of the mesh
        :param dimensions: The absolution dimensions of the mesh
        :return: A constructed finite difference mesh
        """
        self.num_ghost = num_ghost

        xdims = np.linspace(0, dimensions[0], shape[0])
        ydims = np.linspace(0, dimensions[1], shape[1])

        xdiff = np.diff(xdims)
        ydiff = np.diff(ydims)

        for ghost_no in xrange(num_ghost):
            xdims = np.insert(xdims, 0, xdims[0] - xdiff[0])
            xdims = np.append(xdims, xdims[-1] + xdiff[-1])
            ydims = np.insert(ydims, 0, ydims[0] - ydiff[0])
            ydims = np.append(ydims, ydims[-1] + ydiff[-1])

        self.xnodes, self.ynodes = np.meshgrid(xdims, ydims, indexing='ij')

        self.shape = shape
        self.dimensions = dimensions
        self.fields = {}
        self.boundary_fields = {}

    def interior_shape(self):
        return self.shape[0] - 2*self.num_ghost, self.shape[1] - 2*self.num_ghost

    def add_fields(self, *args):
        """
        Add fields to the finite difference mesh
        :param args: The names of the fields to be added
        :return:
        """
        for field_name in args:
            self.fields[field_name] = np.zeros(self.shape, order='F')
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
        nghost = self.num_ghost
        return tuple([self.fields[field_name][nghost:-nghost, nghost:-nghost] for field_name in args])

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
        return self.xnodes[i - self.num_ghost, j - self.num_ghost], self.ynodes[i - self.num_ghost, j - self.num_ghost]

    def nodes(self):
        """
        :return: A tuple containing all of the node coordinates
        """
        return [[(x, y) for x, y in zip(colx, coly)] for colx, coly in zip(self.xnodes, self.ynodes)]

    def plot(self, **kwargs):
        plt.plot(self.xnodes, self.ynodes,
                 np.transpose(self.xnodes), np.transpose(self.ynodes),
                 color='gray', linewidth=1.5)

        if kwargs.get('mark_nodes', False):
            nghost = self.num_ghost
            xnodes, ynodes = self.xnodes[nghost:-nghost, nghost:-nghost], self.ynodes[nghost:-nghost, nghost:-nghost]

            plt.plot(xnodes, ynodes,
                     's', color='black', markersize=3)

            if kwargs.get('show_ghost', False):
                xnodes, ynodes = self.xnodes[0:nghost, 0:], self.ynodes[0:nghost, 0:]
                xnodes, ynodes = np.append(xnodes, self.xnodes[-nghost:, 0:]), np.append(ynodes, self.ynodes[-nghost:, 0:])
                plt.plot(xnodes, ynodes,
                         's', color='blue', markersize=3)
                xnodes, ynodes = self.xnodes[nghost:-nghost, 0:nghost], self.ynodes[nghost:-nghost, 0:nghost]
                xnodes, ynodes = np.append(xnodes, self.xnodes[nghost:-nghost, -nghost:]), np.append(ynodes, self.ynodes[nghost:-nghost, -nghost:])
                plt.plot(xnodes, ynodes,
                         's', color='blue', markersize=3)

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
    grid = FdGrid2D((21, 21), (1, 1), 2)
    grid.add_fields('u', 'v', 'p', 'm', 'pCorr')

    grid.plot(mark_nodes=True, show_ghost=True)
    grid.plot_show()