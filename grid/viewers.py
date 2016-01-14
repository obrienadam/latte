import matplotlib
matplotlib.use('Qt4Agg')

import matplotlib.pylab as plt
import numpy as np
from io import BytesIO


def display_grid(grid, **kwargs):
    fig = plt.figure()
    plt.axes().set_aspect('equal')

    if kwargs.get('mark_core_cells', True):
        core_cell_coords = grid._cell_nodes[1:-1, 1:-1]
        cellx, celly = core_cell_coords[:, :, 0], core_cell_coords[:, :, 1]
        plt.plot(cellx, celly, '-o', np.transpose(cellx), np.transpose(celly), '-o', color='red')

    if kwargs.get('mark_boundary_cells', True):
        boundary_cell_coords = grid._cell_nodes[0, :], \
                               grid._cell_nodes[-1, :], \
                               grid._cell_nodes[1:-1, 0], \
                               grid._cell_nodes[1:-1, -1]

        for coords in boundary_cell_coords:
            plt.plot(coords[:, 0], coords[:, 1], '-x', color='blue')

    if kwargs.get('show', False):
        plt.show()

    f = BytesIO()
    plt.savefig(f)

    return f

def display_fv_solution(grid, field_name, **kwargs):
    fig = plt.figure()
    plt.axes().set_aspect('equal')

    var, = grid.get_cell_fields(field_name)
    cell_nodes = grid.cell_nodes

    plt.contourf(cell_nodes[:, :, 0], cell_nodes[:, :, 1], var, kwargs.get('contour_levels', 15))

    if kwargs.get('show_grid', False):
        plt.plot(cell_nodes[:, :, 0], cell_nodes[:, :, 1], '-', color='black')
        plt.plot(np.transpose(cell_nodes[:, :, 0]), np.transpose(cell_nodes[:, :, 1]), '-', color='black')

    if kwargs.get('mark_cell_centers', False):
        plt.plot(cell_nodes[:, :, 0], cell_nodes[:, :, 1], 'o', color='black')

    if kwargs.get('show', False):
        plt.show()

    f = BytesIO()
    plt.savefig(f)

    return f

def plot_line(grid, field_name, coord, **kwargs):
    fig = plt.figure()
    plt.axes().set_aspect('equal')

    var, = grid.get_cell_fields(field_name)
    cell_nodes = grid.cell_nodes

    axis = kwargs.get('axis', 0)

    if axis == 0:
        plt.plot(cell_nodes[:, coord, 0], var[:, coord])
    elif axis == 1:
        plt.plot(cell_nodes[coord, :, 1], var[coord, :])

    if kwargs.get('show', False):
        plt.show()

    f = BytesIO()
    plt.savefig(f)

    return f

def quiver(grid, u_field_name, v_field_name, **kwargs):
    fig = plt.figure()
    plt.axes().set_aspect('equal')

    u, v = grid.get_cell_fields(u_field_name, v_field_name)
    cell_nodes = grid.cell_nodes

    plt.quiver(cell_nodes[:, :, 0], cell_nodes[:, :, 1], u, v)

    if kwargs.get('show_grid', False):
        plt.plot(cell_nodes[:, :, 0], cell_nodes[:, :, 1], '-', color='black')
        plt.plot(np.transpose(cell_nodes[:, :, 0]), np.transpose(cell_nodes[:, :, 1]), '-', color='black')

    if kwargs.get('mark_cell_centers', False):
        plt.plot(cell_nodes[:, :, 0], cell_nodes[:, :, 1], 'o', color='black')

    if kwargs.get('show', False):
        plt.show()

    f = BytesIO()
    plt.savefig(f)

    return f