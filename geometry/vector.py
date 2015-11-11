def dot(u, v):
    return u[0]*v[0] + u[1]*v[1]

class Vector2D(object):
    def __init__(self, x=0., y=0.):
        self.x, self.y = float(x), float(y)

    def __getitem__(self, item):
        return (self.x, self.y)[item]

