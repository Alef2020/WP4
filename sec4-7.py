import math
from fastener_calculations import *

#def sc_wall_bearing_check(design):
#    failure_list = []
#    for fastener in design.FastenerList:
#        for LC in fastener.load_cases:
#            P_i = math.sqrt(LC.Force_X**2 + LC.Force_Z**2)
#            sigma_lug = P_i / (fastener.configuration.d_s * design.t_2)
#            sigma_wall = P_i / (fastener.configuration.d_s * design.t_3)
#            if sigma_lug > design.material.strength:
#                failure_list.append((fastener.position), "Lug Failure")
#            elif sigma_wall > design.material_wall.strength:
#                failure_list.append((fastener.position), "Wall Failure") 
#    
#    if len(failure_list) > 0:
#        return False
#    else:
#        return True

Total_load_cases = [Load_Case(490.5, 490.5, 1373.4, 137.3 ,24.53 ,49.05)]

pos = [
    [-0.0175, 0, 0.005],
    [0.0175, 0, 0.005],
    [-0.0175, 0, -0.005],
    [0.0175, 0, -0.005]
]

design = design_Configuration(pos, Total_load_cases, 0.05, 0.025, 0.002, 0.002, m2014_T6, Alumunium, Stainless_M5X30)

print(design.FastenerList[3].load_cases[0].Force_Y)