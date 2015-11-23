def compute_area(*args):
    """
    :param args: Points represented as tuples, assumed ordered counter-clockwise
    :return: The area of the polygon as a float
    """
    area = 0.
    num_pts = len(args)

    assert num_pts > 2

    for i in xrange(num_pts):
        area += (args[(i + 1)%num_pts][0] + args[i][0])*(args[(i + 1)%num_pts][1] - args[i][1])

    return 0.5*area

def compute_centroid(*args):
    """
    :param args: Points represented as tuples, assumed ordered counter-clockwise
    :return: The (x, y) coordinates of the centroid
    """
    num_pts = len(args)
    areas = []
    centroids = []

    assert num_pts > 2

    pt0 = args[0]

    for i in xrange(2, num_pts):
        pt1 = args[i - 1]
        pt2 = args[i]
        areas.append(compute_area(pt0, pt1, pt2))
        centroid = (pt0[0] + pt1[0] + pt2[0])/3., (pt0[1] + pt1[1] + pt2[1])/3.
        centroids.append(centroid)

    total_area = sum(areas)
    weights = [area/total_area for area in areas]

    return sum([w*c[0] for w, c in zip(weights, centroids)]), sum([w*c[1] for w, c in zip(weights, centroids)])



if __name__ == '__main__':
    pts = (0.,0.), (1,0.), (1., 1.), (0, 1)

    print 'Area:', compute_area(*pts)
    print 'Centroid: ', compute_centroid(*pts)