import math
from scipy.optimize import minimize
import numpy as np
from fastener_calculations import *

class Bearing_force: # Defines the in plane forces acting on each fastener and calculates the resultant force on the fastener
    def __init__(self, fx, fz):# From what I understand, the forces Fx and Fy from the fastener_calculations already take into accound the portential force from any potential moments
        self.fx = fx 
        self.fz = fz

    def resultant_force(self): # This is the total in plane force on each fastener
        P1 = math.sqrt(self.fx**2 + self.fz**2)
        return P1

# Example of getting the resultant force on each fastener for a desing configuration
design = design_Configuration(Fastener_positions, Total_load_cases, fastener_config = M5_steel)
resultant_forces = [] # position, load case name, resultant force
for fastener in design.FastenerList:
    for LC in fastener.load_cases:
        bf = Bearing_force(LC.Force_X, LC.Force_Z)
        P_i = bf.resultant_force()
        resultant_forces.append((fastener.position, LC.name, P_i))
    

#Does it meet the bearing failure criteria?
class Bearing_failure: # Determines if the fastener will fail or not for the given dimensions
    def __init__(self, P1, D2, t2):
        self.P1 = P1
        self.D2 = D2
        self.t2 = t2

    def bearing_stress(self): # Calculates the bearing stress on the fastener
        sigma_bearing = self.P1 / (self.D2 * self.t2)
        return sigma_bearing
    
    def is_bearing_failure(self, material_bearing_strength): # Checks if the bearing stress exceeds the material bearing strength
        failure_criteria = self.bearing_stress() > self.material_bearing_strength
        if failure_criteria == True:
            return "Failure"
        else:
            return "No Failure"

# Example of checking bearing failure using the same design as before
steel = MaterialProperties(E = 200e9, rho = 7850, strength = 1.67e8) # Example material  # Note that this should be above the design = ...
M5_steel = Fastener_Configuration(0.0083, 0.0075, 0.005, steel) # Example fastener config # Note that this should be above the design = ...
bearing_failure_results = [] # position, failure result
for fastener in design.FastenerList:
    for LC in fastener.load_cases:
        bf = Bearing_force(LC.Force_X, LC.Force_Z)
        P_i = bf.resultant_force()
        bearing_check = Bearing_failure(P_i, fastener.configuration.d_o, design.t_2)
        bearing_check.material_bearing_strength = fastener.configuration.material.strength
        result = bearing_check.is_bearing_failure(fastener.configuration.material.strength)
        bearing_failure_results.append((fastener.position, result))

# Bearing check must also be done for the spacecraft wall 
class SC_wall_properties: # Defines the properties of the spacecraft wall
    def __init__(self, t3, material = Alumunium): 
        self.t3 = t3 # Note that t3 is a local thickness of the spacecraft wall and that the thickness may be different elsewhere
        self.material = material

class SC_wall_bearing_failure(Bearing_failure): # Inherits from Bearing_failure to check bearing failure on the spacecraft wall
    def __init__(self, P1, D2, t3, sc_wall_properties):
        super().__init__(P1, D2, t3)
        self.sc_wall_properties = sc_wall_properties

    def is_sc_wall_bearing_failure(self): # Checks if the bearing stress exceeds the spacecraft wall material bearing strength
        failure_criteria = self.bearing_stress() > self.sc_wall_properties.material.strength
        if failure_criteria == True:
            return "SC Wall Failure"
        else:
            return "SC Wall No Failure"