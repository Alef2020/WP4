import numpy as np
import math
from fastener_calculations import *
from scipy.optimize import minimize

#define loadcases to be analysed
Total_load_cases = [Load_Case(490.5, 490.5, 1373.4, 137.3 ,24.53 ,49.05)] #should eventually be a list of all possible load cases from 4.1

#function for pull through check
def pull_through_check(params):
    design_option = get_design(params)
    for fastener in design_option.FastenerList:
        config = fastener.configuration
        F_1 = math.pi*config.d_h*design_option.t_2 * design_option.material.shear #Force required to pull through head through lug plate
        F_2 = math.pi*config.d_n*design_option.t_3 * design_option.material_wall.shear #Force required to pull through nut through mounting plate
        F_max  = min(F_1, F_2)
        for LC in fastener.load_cases:
            if F_max > abs(LC.Force_Y):
                return 1
    return 0

#This function checks for bearing deformation
def bearing_check(params):
    design = get_design(params)
    failure_list = []
    for fastener in design.FastenerList:
        for LC in fastener.load_cases:
            P_i = math.sqrt(LC.Force_X**2 + LC.Force_Z**2)
            sigma_lug = P_i / (fastener.configuration.d_s * design.t_2)
            sigma_wall = P_i / (fastener.configuration.d_s * design.t_3)
            if sigma_lug > design.material.stress_yield:
                failure_list.append((fastener.position, "Lug Failure"))
            elif sigma_wall > design.material_wall.stress_yield:
                failure_list.append((fastener.position, "Wall Failure")) 

    if len(failure_list) > 0:
        return 1
    else:
        return 0
    
#To get the maximum cost for 
def Max_cost_condition(design):
    if design.cost < 100:
        return 1
    else:
        return 0

#def get_design(params):
#    w ,h, t2 = params 
#    e = 0.3 #distance from edge to center of fastener should be defined based on fastener size
#    fastener_positions = [
#        [ w/2-e, 0,  h/2-e],
#        [-w/2+e, 0,  h/2-e],
#        [ w/2-e, 0, -h/2+e],
#        [-w/2+e, 0, -h/2+e],
#    ]
#    design_option = design_Configuration(fastener_positions, Total_load_cases, width = w, height = h, t_2 = t2, t_3 = 0.002, fastener_config = M5_steel)
#    return design_option
#
##This code unfortunately sucks, I made a mistake in the beginning about how the optimization script works
#def f(params):
#    design_option = get_design(params)
#    return design_option.mass()
#
#constraints = [
#    {'type': 'eq', 'fun': bearing_check},
#    {'type': 'eq', 'fun': pull_through_check}
#]
#
#result = minimize( 
#    f, 
#    [0.5, 0.3, 0.1],
#    bounds=[(0.08 , 10), (0.05,10), (0.002,10)], 
#    method='COBYLA', 
#    constraints=constraints,
#    options={'maxiter': 100, 'disp': True}
#    )
#
#par = result.x[0], result.x[1], result.x[2] 
#print(get_design(par).mass())
#
#print("Optimal design parameters (width, height, t2):", result.x)
#print("Code completed with success:", result.success)
#print("Code message:", result.message)

Results = []
for bolt in Bolt_lst:
    mat_result = []
    for Mat in Material_lst:
        def get_design(params):
            w ,h, t2 = params 
            e = bolt.d_s*1.5 #distance from edge to center of fastener should be defined based on fastener size
            fastener_positions = [
                [ w/2-e, 0,  h/2-e],
                [-w/2+e, 0,  h/2-e],
                [ w/2-e, 0, -h/2+e],
                [-w/2+e, 0, -h/2+e],
            ]
            design_option = design_Configuration(fastener_positions, Total_load_cases, width = w, height = h, t_2 = t2, t_3 = 0.002, Material=Mat, fastener_config = bolt)
            return design_option

        #This code unfortunately sucks, I made a mistake in the beginning about how the optimization script works
        def f(params):
            design_option = get_design(params)
            return design_option.mass()

        constraints = [
            {'type': 'eq', 'fun': bearing_check},
            {'type': 'eq', 'fun': pull_through_check}
        ]

        result = minimize( 
            f, 
            [0.5, bolt.d_s*5.5, 0.1],
            bounds=[(0.03 , 10), (bolt.d_s*5,10), (0.002,10)], 
            method='COBYLA', 
            constraints=constraints,
            options={'disp': True}
            )
        
        par = result.x[0],result.x[1],result.x[2]
        mass = get_design(par).mass()
        mat_result.append([result.x, mass,result.success])
    Results.append(mat_result)


for result_row in Results:
    row = ""
    for optimal in result_row:
        row += f" {(optimal[1]*1000):.2f}, " if optimal[2] else " - "
    print(row)

run = False
while run:
    ft = int(input("Info about what bolt type: "))
    Mt = int(input("what material type: "))

    print(f'Optimal dimensions are: w = {Results[ft][Mt][0][0]}, h = {Results[ft][Mt][0][1]}, t = {Results[ft][Mt][0][2]}')
    print(f'Resulting mass: {Results[ft][Mt][1]} g')
    
    run = True if input("Continue? y/n") == "y" else False