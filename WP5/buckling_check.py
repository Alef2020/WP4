import math
def get_euler_column_buckling_stress(youngs_modulus,area_moment_of_inertia,area,length):
    return math.pi**2 * youngs_modulus * area_moment_of_inertia/(area*length**2)

def get_euler_column_buckling_safety_margin(youngs_modulus,area_moment_of_inertia,area,length,applied_force,):
    applied_stress = applied_force/area
    return applied_stress/get_euler_column_buckling_stress(youngs_modulus,area_moment_of_inertia,area,length)

def get_shell_buckling_stress(youngs_modulus,length,pressure,radius,thickness,internal_pressure,poisson_ratio)
    k = _lambda + 12/(math.pi**4) * length**4/(radius**2 * thickness**2) * (1-poisson_ratio**2)/_lambda
