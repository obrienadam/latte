import numpy as np
from finite_difference import FdGrid2D
import matplotlib.pyplot as plt
from geometry import polygon as poly


class FvGrid2D(FdGrid2D):
    def __init__(self, shape, dimensions, num_ghost=1):
        """
        Initialize a 2D finite volume grid
        :param shape: The nodal resolution of the mesh
        :param dimensions: The absolute dimensions of the mesh
        :return: A constructed finite volume mesh
        """
        super(FvGrid2D, self).__init__(shape, dimensions, num_ghost)

        xdims = np.linspace(0, dimensions[0], shape[0])
        ydims = np.linspace(0, dimensions[1], shape[1])

        xdiff = np.diff(xdims)
        ydiff = np.diff(ydims)

        for ghost_no in xrange(num_ghost):
            xdims = np.insert(xdims, 0, xdims[0] - xdiff[0])
            xdims = np.append(xdims, xdims[-1] + xdiff[-1])
            ydims = np.insert(ydims, 0, ydims[0] - ydiff[0])
            ydims = np.append(ydims, ydims[-1] + ydiff[-1])

        deltaxs = np.diff(xdims)
        deltays = np.diff(ydims)

        xcdims = [xdim + 0.5 * deltax for xdim, deltax in zip(xdims, deltaxs)]
        ycdims = [ydim + 0.5 * deltay for ydim, deltay in zip(ydims, deltays)]

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
        self.hface_norms_x = -np.diff(self.ynodes, axis=0)
        self.hface_norms_y = np.diff(self.xnodes, axis=0)

        self.vface_norms_x = np.diff(self.ynodes, axis=1)
        self.vface_norms_y = -np.diff(self.xnodes, axis=1)

        # Compute diffusion metrics

        self.cell_shape = len(xcdims), len(ycdims)
        self.hface_shape = len(xcdims), len(ydims)
        self.vface_shape = len(xdims), len(ycdims)

    def interior_shape(self):
        return self.cell_shape[0] - 2*self.num_ghost, self.cell_shape[1] - 2*self.num_ghost

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

            field['cell_data'] = np.zeros(self.cell_shape, order='F') if kwargs.get('cell_centered', False) else None
            field['node_data'] = np.zeros(self.shape, order='F') if kwargs.get('node_centered', False) else None
            field['hface_data'] = np.zeros(self.hface_shape, order='F') if kwargs.get('hface_centered', False) else None
            field['vface_data'] = np.zeros(self.vface_shape, order='F') if kwargs.get('vface_centered', False) else None

            if field['cell_data'] != None:
                cell_data = field['cell_data']
                ng = self.num_ghost
                field['boundary_data'] = {'east': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[-ng:, :]},
                                          'west': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[0:ng, :]},
                                          'north': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[:, -ng:]},
                                          'south': {'type': 'fixed', 'refval': 0., 'bfield': cell_data[:, 0:ng]}}

    def rename_field(self, **kwargs):
        """
        Renames the specifies field
        :param kwargs: Key-value pair where the key is the old name and the value is the new name
        :return:
        """
        for old_field_name in kwargs:
            self.fields[kwargs[old_field_name]] = self.fields.pop(old_field_name)

    def set_const_field(self, **kwargs):
        """
        Sets the specified fields to a constant value
        :param kwargs: Field names and their values
        :return:
        """
        for field_name in kwargs:
            field = self.fields[field_name]

            if not field['cell_data'] == None:
                field['cell_data'][:, :] = kwargs[field_name]

            if not field['node_data'] == None:
                field['node_ata'][:, :] = kwargs[field_name]

            if not field['hface_data'] == None:
                field['hface_data'][:, :] = kwargs[field_name]

            if not field['vface_data'] == None:
                field['vface_data'][:, :] = kwargs[field_name]


    def cell_data(self, *args):
        """
        :param args: Names of fields
        :return: A tuple containing the cell data of the desired fields
        """
        ng = self.num_ghost
        return tuple([self.fields[field_name]['cell_data'][ng:-ng, ng:-ng] for field_name in args])

    def boundary_data(self, *args):
        """
        :param args: Name of fields
        :return: A tuple containing the cell data of the desired fields
        """
        bcs = []
        for field_name in args:
            bcs.append(self.fields[field_name]['boundary_data'])

        return tuple(bcs)

    def node_data(self, *args):
        """
        :param args: Names of fields
        :return: A tuple containing the node data of the desired fields
        """
        nghost = self.num_ghost
        return tuple([self.fields[field_name]['node_data'][nghost:-nghost, nghost:-nghost] for field_name in args])

    def face_data(self, field_name):
        """
        :param args: Names of fields
        :return: A tuple containing tuples with the horizontal and vertical face data of the desired fields
        """
        ng = self.num_ghost
        hf, vf = self.fields[field_name]['hface_data'], self.fields[field_name]['vface_data']
        return vf[ng + 1:-ng, ng:-ng], vf[ng:-ng - 1, ng:-ng], hf[ng:-ng, ng + 1:-ng], hf[ng:-ng, ng:-ng - 1]

    def xc(self, i, j):
        """
        :param i: Index i
        :param j: Index j
        :return: Tuple containing the x, y coordinate of the cell center
        """
        return self.cell_centers_x[i, j], self.cell_centers_y[i, j]

    def face_norms(self):
        return self.east_face_norms(), self.west_face_norms(), self.north_face_norms(), self.south_face_norms()

    def east_face_norms(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.vface_norms_x[1 + ng:-ng, ng:-ng],
                                                                   self.vface_norms_y[1 + ng:-ng, ng:-ng])])

    def west_face_norms(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(-self.vface_norms_x[ng:-ng - 1, ng:-ng],
                                                                   self.vface_norms_y[ng:-ng - 1, ng:-ng])])

    def north_face_norms(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.hface_norms_x[ng:-ng, 1 + ng:-ng],
                                                                   self.hface_norms_y[ng:-ng, 1 + ng:-ng])])

    def south_face_norms(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.hface_norms_x[ng:-ng, ng:-ng - 1],
                                                                   -self.hface_norms_y[ng:-ng, ng:-ng - 1])])

    def cell_rvecs(self):
        return self.east_cell_rvecs(), self.west_cell_rvecs(), self.north_cell_rvecs(), self.south_cell_rvecs()

    def east_cell_rvecs(self):
        ng = self.num_ghost
        if ng > 1:
            return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.rcell_h_x[ng:-ng + 1, ng:-ng],
                                                                       self.rcell_h_y[ng:-ng + 1, ng:-ng])])
        else:
            return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.rcell_h_x[ng:, ng:-ng],
                                                                       self.rcell_h_y[ng:, ng:-ng])])

    def west_cell_rvecs(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(-self.rcell_h_x[ng - 1:-ng, ng:-ng],
                                                                   -self.rcell_h_y[ng - 1:-ng, ng:-ng])])

    def north_cell_rvecs(self):
        ng = self.num_ghost
        if ng > 1:
            return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.rcell_v_x[ng:-ng, ng:-ng + 1],
                                                                       self.rcell_v_y[ng:-ng, ng:-ng + 1])])
        else:
            return np.array([zip(xvals, yvals) for xvals, yvals in zip(self.rcell_v_x[ng:-ng, ng:],
                                                                       self.rcell_v_y[ng:-ng, ng:])])

    def south_cell_rvecs(self):
        ng = self.num_ghost
        return np.array([zip(xvals, yvals) for xvals, yvals in zip(-self.rcell_v_x[ng:-ng, ng - 1:-ng],
                                                                   -self.rcell_v_y[ng:-ng, ng - 1:-ng])])

    def cell_centers(self):
        """
        :return: A tuple containing all the coordinates of the cell centers
        """
        ng = self.num_ghost
        return self.cell_centers_x[ng:-ng, ng:-ng], self.cell_centers_y[ng:-ng, ng:-ng]

    def plot(self, **kwargs):
        super(FvGrid2D, self).plot(**kwargs)

        nghost = self.num_ghost

        if kwargs.get('mark_cell_centers', False):
            if kwargs.get('show_ghost', False):
                cell_centers_x, cell_centers_y = self.cell_centers_x, self.cell_centers_y
            else:
                cell_centers_x, cell_centers_y = self.cell_centers_x[nghost:-nghost, nghost:-nghost], \
                                                 self.cell_centers_y[nghost:-nghost, nghost:-nghost]

            plt.plot(cell_centers_x, cell_centers_y,
                     'o', color='red', markersize=3)

        if kwargs.get('mark_faces', False):
            if kwargs.get('show_ghost', False):
                hface_centers_x, hface_centers_y = self.hface_centers_x, self.hface_centers_y
                vface_centers_x, vface_centers_y = self.vface_centers_x, self.vface_centers_y
            else:
                hface_centers_x, hface_centers_y = self.hface_centers_x[nghost:-nghost, nghost:-nghost], \
                                                   self.hface_centers_y[nghost:-nghost, nghost:-nghost]
                vface_centers_x, vface_centers_y = self.vface_centers_x[nghost:-nghost, nghost:-nghost], \
                                                   self.vface_centers_y[nghost:-nghost, nghost:-nghost]

            plt.plot(hface_centers_x, hface_centers_y,
                     'p', color='blue', markersize=3)
            plt.plot(vface_centers_x, vface_centers_y,
                     'p', color='blue', markersize=3)


if __name__ == '__main__':
    grid = FvGrid2D((41, 41), (1, 1), 1)
    grid.add_fields('u', 'v', 'p', 'pCorr', 'm', 'phi', cell_centered=True, face_centered=True)

    sfe = grid.east_face_norms()
    sfw = grid.west_face_norms()

    grid.plot(mark_nodes=True, mark_cell_centers=True, mark_faces=True)
    grid.plot_show()
