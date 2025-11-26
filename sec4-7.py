import math
from fastener_calculations import *

def sc_wall_bearing_check(design):
    failure_list = []
    for fastener in design.FastenerList:
        for LC in fastener.load_cases:
            P_i = math.sqrt(LC.Force_X**2 + LC.Force_Z**2)
            sigma_lug = P_i / (fastener.configuration.d_s * design.t_2)
            sigma_wall = P_i / (fastener.configuration.d_s * design.t_3)
            if sigma_lug > design.material.strength:
                failure_list.append((fastener.position), "Lug Failure")
            elif sigma_wall > design.material_wall.strength:
                failure_list.append((fastener.position), "Wall Failure") 
    
    if len(failure_list) > 0:
        return False
    else:
        return True