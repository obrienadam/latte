from boundary_condition import BoundaryCondition

class FixedValue(BoundaryCondition):
    def __init__(self, ref, **kwargs):
        self.ref = ref
        self.value = kwargs['value']

    def set(self):
        self.ref[:, :] = self.value
