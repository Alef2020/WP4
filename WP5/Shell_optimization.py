#import required libraries and functions
import math
from materials import * 
from buckling_check import *
from scipy.optimize import minimize

def f(params, rho = material.density):
    R, t, L = params
    volume = 2*R*math.pi*t*L
    mass = volume * rho
    return mass

conditions = [
    #{'type':'ineq', 'fun': get_euler_column_buckling_safety_margin},
    #{'type':'ineq', 'fun': get_shell_buckling_safety_margin},
    {'type':'ineq', 'fun': volume_constraint},
    {'type':'ineq', 'fun': tank_integration_constraint}
]

result = minimize(
f, 
[0.597, 0.002, 5], 
constraints = conditions,
bounds=[(0.2, 1),(0.001,0.02),(3, 7)], 
method='SLSQP',
options={'disp':True}
)

print(result.success)
print("the optimal design parameters are:")
print(f'radius: {result.x[0]} m')
print(f'radius: {result.x[1]*1000} mm')
print(f'radius: {result.x[2]} m')