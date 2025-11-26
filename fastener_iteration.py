import numpy as np
import math
from fastener_calculations import *
from scipy.optimize import minimize

#define loadcases to be analysed
Total_load_cases = [Load_Case(3,4,2,3,6,1), Load_Case(9,2,4,1,3,4)] #should eventually be a list of all possible load cases from 4.1


#function of checking the spacing constraint of fasteners
def check_spacing_constraint(design_option):#checking constraints for the spacing of the fasteners
    #positions = np.array(fastener_positions, dtype=float)
    #min_dist = 2 * D2
    #max_dist = 3 * D2
    #n = len(positions)
    #for i in range(n):
    #    for j in range(i+1, n):
    #        
    #        dz = abs(positions[i][2] - positions[j][2])
    #        if  dz < min_dist  or dz > max_dist:
    #            return 0
    return 0

#function for pull through check
def pull_through_check(design_option):
    for fastener in design_option.FastenerList:
        config = fastener.configuration
        F_1 = math.pi * (config.d_h/2)**2 * config.material.shear #Force required to pull through head through lug plate
        F_2 = math.pi * (config.d_n/2)**2 * config.material.shear #Force required to pull through nut through mounting plate
        F_max  = min(F_1, F_2)
        for LC in fastener.load_cases:
            if F_max > LC.Force_Y:
                return 1
    return 0

def bearing_check(design):
    failure_list = []
    for fastener in design.FastenerList:
        for LC in fastener.load_cases:
            P_i = math.sqrt(LC.Force_X**2 + LC.Force_Z**2)
            sigma_lug = P_i / (fastener.configuration.d_s * design.t_2)
            sigma_wall = P_i / (fastener.configuration.d_s * design.t_3)
            if sigma_lug > design.material.strength:
                failure_list.append((fastener.position, "Lug Failure"))
            elif sigma_wall > design.material_wall.strength:
                failure_list.append((fastener.position, "Wall Failure")) 
    
    if len(failure_list) > 0:
        return 1
    else:
        return 0

def get_design(params):
    w ,h, t2, t3 = params 
    e = 0.3 #distance from edge to center of fastener should be defined based on fastener size
    fastener_positions = [
        [ w/2-e, 0,  h/2-e],
        [-w/2+e, 0,  h/2-e],
        [ w/2-e, 0, -h/2+e],
        [-w/2+e, 0, -h/2+e],
    ]
    design_option = design_Configuration(fastener_positions, Total_load_cases, width = w, height = h, t_2 = t2, t_3 = t3, fastener_config = M5_steel)
    return design_option


def Mass_function(params):
    w ,h, t2, t3 = params 
    e = 0.3 #distance from edge to center of fastener should be defined based on fastener size
    fastener_positions = [
        [ w/2-e, 0,  h/2-e],
        [-w/2+e, 0,  h/2-e],
        [ w/2-e, 0, -h/2+e],
        [-w/2+e, 0, -h/2+e],
    ]

    design_option = design_Configuration(fastener_positions, Total_load_cases, width = w, height = h, t_2 = t2, t_3 = t3, fastener_config = M5_steel)
    return design_option.mass()


constraints = [
    {'type': 'eq', 'Mass_function': lambda params: pull_through_check(get_design(params))}
    #{'type': 'eq', 'Mass_function': lambda params: bearing_check(get_design(params))},
    #{'type': 'eq', 'Mass_function': lambda params: check_spacing_constraint(get_design(params))}
]

result = minimize(
    Mass_function, 
    [3, 2, 4, 6], 
    bounds=[(2,5), (1,4), (2,6), (4,8)], 
    method='SLSQP', 
    constraints=constraints,
    options={'maxiter': 20}
    )

print("Optimal design parameters (width, height, t2, t3):", result.x, "   m:  ", result.Mass_function)
print("Code completed with success:", result.success)