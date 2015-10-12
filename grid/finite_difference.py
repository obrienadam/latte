import numpy as np

class FdGrid2D(object):
    def __init__(self, shape, dimensions):
        xdims = np.linspace(0, dimensions[0], shape[0])
        ydims = np.linspace(0, dimensions[1], shape[1])
        self.xnodes, self.ynodes = np.meshgrid(xdims, ydims, indexing='ij')

        self.shape = shape
        self.fields = {}
        self.boundary_fields = {}

    def add_fields(self, *args):
        for field_name in args:
            self.fields[field_name] = np.zeros(self.shape)
            field = self.fields[field_name]
            self.boundary_fields[field_name] = field[-1, :], field[0, :], field[:, -1], field[:, 0]

    def get_field_names(self):
        return self.fields.keys()

    def node(self, index):
        return self.xnodes[index], self.ynodes[index]

    def nodes(self):
        return self.xnodes, self.ynodes

if __name__ == '__main__':
    grid = FdGrid2D((41, 21), (1, 1))
    grid.add_fields('u', 'v', 'p', 'm', 'pCorr')

    print grid.node((40, 20))