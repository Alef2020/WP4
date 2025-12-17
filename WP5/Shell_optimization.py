#import required libraries and functions
import math
from materials import * 
from buckling_check import *
from scipy.optimize import minimize



cparams=0.597, 0.002, 5
check = True
message = ""
if get_euler_column_buckling_safety_margin(cparams) < 0:
    pass
    check = False
    #message += "euler_buckling"
elif get_shell_buckling_safety_margin(cparams) <0:
    check = False
    message += "shell_buckling"
elif volume_constraint(cparams) <0:
    check = False
    message += "volume"
elif tank_integration_constraint(cparams) <0:
    check = False
    message += "integration"

print("The starting values comply with the constraints:")
if not check:print(message)

def f(params, rho = material.density):
    R, t, L = params
    volume = 2*R*math.pi*t*L
    mass = volume * rho
    return mass

conditions = [
    #{'type':'ineq', 'fun': get_euler_column_buckling_safety_margin},
    {'type':'ineq', 'fun': get_shell_buckling_safety_margin},
    {'type':'ineq', 'fun': volume_constraint},
    {'type':'ineq', 'fun': tank_integration_constraint}
]

result = minimize(
f, 
[0.597, 0.002, 5], 
constraints = conditions,
bounds=[(0.01, None),(0.001,None),(0.01,None)], 
method='SLSQP',
#options={'disp':True}
)

#print(result.success)
#print("the optimal design parameters are:")
#print(f'radius: {result.x[0]} m')
#print(f'thickness: {result.x[1]*1000} mm')
#print(f'length: {result.x[2]} m')
#
#print(result.message)