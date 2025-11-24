import math
from scipy.optimize import minimize
import numpy as np

#Resultant Force

class bearingforce:
    def __init__(self, fx, fy, fz):
        self.fx = fx
        self.fy = fy
        self.fz = fz

    def resultant_force(self):
        P1 = math.sqrt(self.fx**2 + self.fy**2 + self.fz**2)
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

