# Buckling check
# By: 6126324
#------------------------


#importing packages
import math
import scipy.optimize as opt
import colorama
colorama.init()
global area_moment_of_inertia




#defining functions as specified in WP5
def get_euler_column_buckling_stress(youngs_modulus,area_moment_of_inertia,area,length):
    return math.pi**2 * youngs_modulus * area_moment_of_inertia/(area*length**2)

def get_shell_buckling_stress(youngs_modulus,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio):
    def get_k(_lambda):
        return _lambda + 12/(math.pi**4) * length**4/(radius**2 * thickness**2) * (1-poisson_ratio**2)/_lambda
    k = opt.minimize_scalar(get_k).fun
    minimum_pressure_difference = internal_pressure - outside_pressure
    q = minimum_pressure_difference/youngs_modulus * (radius/thickness)**2
    return (1.983-0.983*math.e**(-23.14*q))*k*math.pi**2*youngs_modulus/(12*(1-poisson_ratio**2))*(thickness/length)**2

# finding safety margin for both types of buckling
def get_euler_column_buckling_safety_margin(youngs_modulus,area_moment_of_inertia,area,length,applied_force,):
    def check()
    applied_stress = applied_force/area
    return applied_stress/get_euler_column_buckling_stress(youngs_modulus,area_moment_of_inertia,area,length)-1

def get_shell_buckling_safety_margin(youngs_modulus,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio,applied_force):
    applied_stress = applied_force/area
    return applied_stress/get_shell_buckling_stress(youngs_modulus,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio)-1

# gives all buckling information at once in a dictionary
def get_minimum_buckling_safety_margin(youngs_modulus,area_moment_of_inertia,area,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio,applied_force):
    euler_column_buckling_stress = get_euler_column_buckling_stress(youngs_modulus,area_moment_of_inertia,area,length)
    euler_column_buckling_safety_margin = get_euler_column_buckling_safety_margin(youngs_modulus,area_moment_of_inertia,area,length,applied_force)

    shell_buckling_stress = get_shell_buckling_stress(youngs_modulus,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio)
    shell_buckling_safety_margin = get_shell_buckling_safety_margin(youngs_modulus,length,outside_pressure,radius,thickness,internal_pressure,poisson_ratio,applied_force)

    minimum_buckling_stress = min(euler_column_buckling_stress,shell_buckling_stress)
    minimum_buckling_safety_margin = min(euler_column_buckling_safety_margin,shell_buckling_safety_margin)
    buckling_dictionary = {
        "euler column buckling stress": euler_column_buckling_stress,
        "shell buckling stress": shell_buckling_stress,
        "minimum buckling stress": minimum_buckling_stress,
        "euler column buckling safety margin": euler_column_buckling_safety_margin,
        "shell buckling safety margin": shell_buckling_safety_margin,
        "minimum buckling safety margin": minimum_buckling_safety_margin
    }

    #prints minimum safety margin in color
    if minimum_buckling_safety_margin < 0:
        print(Fore.red + "Minimum safety margin: " + str(minimum_buckling_safety_margin) + Style.RESET_ALL)
    else:
        print(Fore.green + "Minimum safety margin " + str(minimum_buckling_safety_margin) + Style.RESET_ALL)
    return buckling_dictionary