import math
from scipy.optimize import minimize
import numpy as np

#Resultant Force

class FastenerLoad:
    def __init__(self, Fx, Fz, FM):
        self.Fx = Fx
        self.Fz = Fz
        self.FM = FM

    def total_inplane_load(self):
        P1 = math.sqrt(self.Fx**2 + self.Fz**2 + self.FM**2)
        return P1

#Does it meet the bearing failure criteria?
class bearingfailure:
    def __init__(self, P1, D2, t2):
        self.P1 = P1
        self.D2 = D2
        self.t2 = t2

    def bearing_stress(self):
        sigma_bearing = self.P1 / (self.D2 * self.t2)
        return sigma_bearing
    def is_bearing_failure(self, material_yield_strength):
        failure = self.bearing_stress() > material_yield_strength
        return failure

