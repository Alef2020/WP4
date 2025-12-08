import math
from fastener_calculations import *


# ===============================
# Assumed material and geometry
# ===============================
thermal_coefficient_lug = 2e-5  # 1/K
thermal_coefficient_bolt = 1e-5  # 1/K
fastener_youngs_modulus = 193e9  # Pa
backplate_youngs_modulus = 50e9  # Pa
skin_youngs_modulus = 50e9  # Pa

backplate_thickness = 0.01  # m
skin_thickness = 0.01  # m
bolt_head_radius = 0.01  # m
hole_diameter = 0.01  # m
bolt_distance_to_edge = 0.2  # m

total_thickness = backplate_thickness + skin_thickness




# ===========================
# Fastener Geometry
# ===========================


d_nominal = 0.005
d_minor = 4.134e-3 # https://www.accu.co.uk/p/117-iso-metric-thread-dimensions

L_h_sub   = 3.5e-3 #m
L_eng_sub  = 3e-3  #m
L_n_sub   =  3.5e-3  #m
L_shank = 4e-3  #m (t2+t3)

d_h_sub = 8e-3 #m
d_eng_sub = 5e-3 #m
d_n_sub = 8e-3 #m



fastener_geometry_list = [
    [L_shank, d_nominal/2],
    [L_h_sub, d_h_sub/2],
    [L_eng_sub, d_eng_sub/2],
    [L_n_sub, d_n_sub/2]
]

integration_steps = 100000  # numerical integration steps

# ===============================
# Temperature conditions
# ===============================
assembly_temperature = 288  # K
min_temp = 263  # K
max_temp = 298  # K

delta_t_list = []
delta_t_list.append(min_temp - assembly_temperature)
delta_t_list.append(max_temp - assembly_temperature)

# ===============================
# Compliance of fastener
# ===============================
compliance_fastener = 0
for fastener in fastener_geometry_list:
    length = fastener[0]
    radius = fastener[1]
    area = math.pi * radius ** 2
    compliance_fastener += length / (fastener_youngs_modulus * area)

# ===============================
# Compliance of clamped part (ESA method)
# ===============================
# Compression cone half-angle (tan θ)
compression_cone_half_angle = 0.362
compression_cone_half_angle += 0.032 * math.log(total_thickness / bolt_head_radius)
compression_cone_half_angle += 0.153 * math.log(bolt_distance_to_edge / (2 * bolt_head_radius))

w = 1  # 1 for nut-tightened joints, 2 for threaded-hole joints

# Limit radius of compression cone
compression_cone_limit_radius = bolt_head_radius + 0.5 * w * total_thickness * compression_cone_half_angle
if compression_cone_limit_radius > bolt_distance_to_edge:
    compression_cone_limit_radius = bolt_distance_to_edge

# Numerical integration for clamped part compliance
step_size = total_thickness / integration_steps
compliance_clamped_part = 0

for i in range(integration_steps):
    z = i * step_size

    # Determine radius at this z
    mid_thickness = total_thickness / 2
    if z <= mid_thickness:
        r = bolt_head_radius + (2 * (compression_cone_limit_radius - bolt_head_radius) / total_thickness) * z
    else:
        r = 2 * compression_cone_limit_radius - bolt_head_radius + (
                    2 * (bolt_head_radius - compression_cone_limit_radius) / total_thickness) * z

    if r > bolt_distance_to_edge:
        r = bolt_distance_to_edge

    # Determine Young's modulus at this z
    if z <= backplate_thickness:
        E = backplate_youngs_modulus
    else:
        E = skin_youngs_modulus

    # Increment compliance
    compliance_clamped_part += step_size / (E * math.pi * r ** 2)

# ===============================
# Calculate forces on fasteners
# ===============================
total_bolt_length = 0
for fastener in fastener_geometry_list:
    total_bolt_length += fastener[0]

for delta_t in delta_t_list:
    force = total_bolt_length * (thermal_coefficient_lug - thermal_coefficient_bolt) * delta_t
    force = force / (compliance_fastener + compliance_clamped_part)
    print("ΔT = ", delta_t, "K, Force = ", force, "N")



# ===============================
# Force Ratio Function
# ===============================

def force_ratio(compliance_fastener, compliance_clamped_part):
    return compliance_fastener / (compliance_fastener + compliance_clamped_part)