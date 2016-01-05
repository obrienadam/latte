import numpy as np

def polygon_area(x_coords, y_coords):
    """
    Computes the area of an arbitrary polygon
    :param x_coords: X-coordinates of polygon, arranged ccw
    :param y_coords: Y-coordinates of polygon, arranged ccw
    :return: Area of the polygon
    """
    npts = len(x_coords)
    area = 0.

    for i in xrange(npts):
        area += x_coords[i]*y_coords[(i + 1)%npts] - x_coords[(i + 1)%npts]*y_coords[i]

    return 0.5*area

def polygon_centroid(x_coords, y_coords):
    """
    :param x_coords: X-coordinates of polygon, arranged ccw
    :param y_coords: Y-coordinates of polygon, arranged ccw
    :return: Numpy array representing the centroid of the polygon
    """
    npts = len(x_coords)
    coeff = 1./(6.*polygon_area(x_coords, y_coords))
    cx, cy = 0., 0.

    for i in xrange(npts):
        ai = x_coords[i]*y_coords[(i + 1)%npts] - x_coords[(i + 1)%npts]*y_coords[i]
        cx += (x_coords[i] + x_coords[(i + 1)%npts])*ai
        cy += (y_coords[i] + y_coords[(i + 1)%npts])*ai

    return coeff*np.array((cx, cy), dtype=float, order='F')

def polygon_face_normals(x_coords, y_coords, start_point=0, nfaces=0):
    """
    :param x_coords: X-coordinates of polygon, arranged ccw
    :param y_coords: Y-coordinates of polygon, arranged ccw
    :param start_point: Set the starting point, numbered ccw
    :param nfaces: Set the number of face normals to be computed. If 0, find all of them
    :return: Numpy array representing outward pointing face normals of the polygon
    """
    npts = len(x_coords)
    nfaces = npts if nfaces is 0 else nfaces
    face_norms = np.zeros((nfaces, 2), dtype=float)

    for i in xrange(start_point, start_point + nfaces):
        dx = x_coords[(i + 1)%npts] - x_coords[i%npts]
        dy = y_coords[(i + 1)%npts] - y_coords[i%npts]

        face_norms[i - start_point] = face_norm(dx, dy)

    return face_norms

def face_norm(dx, dy):
    return dy, -dx

if __name__ == '__main__':
    poly = [0, 1.2, 1, 0], [0, 0, 1, 1]

    print "Area:", polygon_area(*poly)
    print "Centroid:", polygon_centroid(*poly)
    print "Face norms:", polygon_face_normals(*poly)