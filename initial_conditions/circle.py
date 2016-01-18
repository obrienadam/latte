def init_circle(grid, field_name, center, radius, val):
    field, = grid.get_cell_fields(field_name)
    cell_nodes = grid.cell_nodes
    corner_nodes = grid.corner_nodes

    numx, numy = field.shape
    for j in xrange(numy):
        for i in xrange(numx):
            if (cell_nodes[i, j, 0] - center[0])**2 + (cell_nodes[i, j, 1] - center[1])**2 < radius**2:
                field[i, j] = val
            else:
                field[i, j] = 0.