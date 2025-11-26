import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize

material_properties = {
    "2014-T6": {
        "density": 2800,
        "E": 72.4 * (10**9),
        "sigma_y": 414 * (10**6),
        "sigma_u": 483 * (10**6),
    },
    "7075-T6": {
        "density": 2810,
        "E": 71.7 * (10**9),
        "sigma_y": 503 * (10**6),
        "sigma_u": 572 * (10**6),
    },
    "2024-T2": {
        "density": 2780,
        "E": 73.1 * (10**9),
        "sigma_y": 324 * (10**6),
        "sigma_u": 469 * (10**6),
    },
    "2024-T3": {
        "density": 2780,
        "E": 73.1 * (10**9),
        "sigma_y": 345 * (10**6),
        "sigma_u": 483 * (10**6),
    },
    "2024-T4": {
        "density": 2780,
        "E": 71.0 * (10**9),
        "sigma_y": 324 * (10**6),
        "sigma_u": 469 * (10**6),
    },
    "AZ1916-T6": {
        "density": 1810,
        "E": 44.8 * (10**9),
        "sigma_y": 145 * (10**6),
        "sigma_u": 275 * (10**6),
    },
    "356-T6": {
        "density": 2680,
        "E": 72.4 * (10**9),
        "sigma_y": 138 * (10**6),
        "sigma_u": 207 * (10**6),
    },
    "4130-steel": {
        "density": 7850,
        "E": 205 * (10**9),
        "sigma_y": 435 * (10**6),
        "sigma_u": 670 * (10**6),
    },
    "8630-steel": {
        "density": 7850,
        "E": 200 * (10**9),
        "sigma_y": 550 * (10**6),
        "sigma_u": 620 * (10**6),
    },
}


def interpolation(filename, i=0, j=1):
    x_vals = []
    y_vals = []

    with open(filename, "r") as f:
        # read all lines except header
        lines = f.readlines()[1:]

    for line in lines:
        # replace decimal comma and semicolon
        line = line.replace(",", ".").replace(";", " ")
        parts = line.strip().split()

        # skip lines that do not contain enough columns
        if len(parts) <= max(i, j):
            continue

        try:
            x = float(parts[i])
            y = float(parts[j])
            x_vals.append(x)
            y_vals.append(y)
        except ValueError:
            # skip lines with invalid numeric entries
            continue

    # convert to numpy arrays
    x = np.array(x_vals)
    y = np.array(y_vals)

    # create interpolation function
    f = interp1d(x, y, kind="cubic", fill_value="extrapolate")

    return f


def Kt(curve, W_D):
    if curve == "Curve 1":
        f_interp = interpolation("Kt Curves.csv", 0, 1)
    if curve == "Curve 2":
        f_interp = interpolation("Kt Curves.csv", 2, 3)
    if curve == "Curve 3":
        f_interp = interpolation("Kt Curves.csv", 4, 5)
    if curve == "Curve 4":
        f_interp = interpolation("Kt Curves.csv", 6, 7)
    if curve == "Curve 5":
        f_interp = interpolation("Kt Curves.csv", 8, 9)
    if curve == "Curve 6":
        f_interp = interpolation("Kt Curves.csv", 10, 11)
    if curve == "Curve 7":
        f_interp = interpolation("Kt Curves.csv", 12, 13)

    K_t = f_interp(W_D)

    return K_t


def K_Bry(t_D, e_D):
    if t_D <= 0.07:
        f_interp = interpolation("K_bry Curves.csv", 0, 1)
        K_bry = f_interp(e_D)
    elif t_D <= 0.09:
        f_interp = interpolation("K_bry Curves.csv", 2, 3)
        K_bry = f_interp(e_D)
    elif t_D <= 0.11:
        f_interp = interpolation("K_bry Curves.csv", 4, 5)
        K_bry = f_interp(e_D)
    elif t_D <= 0.135:
        f_interp = interpolation("K_bry Curves.csv", 6, 7)
        K_bry = f_interp(e_D)
    elif t_D <= 0.175:
        f_interp = interpolation("K_bry Curves.csv", 8, 9)
        K_bry = f_interp(e_D)
    elif t_D <= 0.25:
        f_interp = interpolation("K_bry Curves.csv", 10, 11)
        K_bry = f_interp(e_D)
    elif t_D <= 0.35:
        f_interp = interpolation("K_bry Curves.csv", 12, 13)
        K_bry = f_interp(e_D)
    elif t_D <= 0.5:
        f_interp = interpolation("K_bry Curves.csv", 14, 15)
        K_bry = f_interp(e_D)
    else:
        f_interp = interpolation("K_bry Curves.csv", 16, 17)
        K_bry = f_interp(e_D)

    return K_bry


def K_ty(curve, Av_Abr):
    if curve == "Curve 1":
        f_interp = interpolation("K_ty Curves.csv", 0, 1)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 2":
        f_interp = interpolation("K_ty Curves.csv", 2, 3)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 3":
        f_interp = interpolation("K_ty Curves.csv", 4, 5)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 4":
        f_interp = interpolation("K_ty Curves.csv", 6, 7)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 5":
        f_interp = interpolation("K_ty Curves.csv", 8, 9)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 6":
        f_interp = interpolation("K_ty Curves.csv", 10, 11)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 7":
        f_interp = interpolation("K_ty Curves.csv", 12, 13)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 8":
        f_interp = interpolation("K_ty Curves.csv", 14, 15)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 9":
        f_interp = interpolation("K_ty Curves.csv", 16, 17)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 10":
        f_interp = interpolation("K_ty Curves.csv", 18, 19)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 11":
        f_interp = interpolation("K_ty Curves.csv", 20, 21)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 12":
        f_interp = interpolation("K_ty Curves.csv", 22, 23)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 13":
        f_interp = interpolation("K_ty Curves.csv", 24, 25)
        Kty = f_interp(Av_Abr)
    if curve == "Curve 14":
        f_interp = interpolation("K_ty Curves.csv", 26, 27)
        Kty = f_interp(Av_Abr)
    return Kty


def load_check(
    material,
    curve_Kt,
    flanges,
    c,
    F_xx,
    F_yy,
    F_zz,
    M_xx,
    M_yy,
    M_zz,
    W,
    D_1,
    t_1,
    e,
    l,
    curve_Kty="Curve 3",
):
    # Av = 6 / (
    #     3 / (0.5 * t_1 * (W - D_1 * np.sin(np.pi / 4)))
    #     + 1 / (0.5 * t_1 * (W - D_1 * np.sin(np.pi / 4)))
    #     + 1 / (0.5 * t_1 * (W - D_1))
    #     + 1 * ((e - D_1 / 2) * t_1)
    # )
    A_1 = t_1 * (W / 2 - D_1 / 2 * np.sin(np.pi / 4))
    A_2 = t_1 * (W / 2 - D_1 / 2)
    A_3 = A_2
    A_4 = A_1
    Av = 6 / ((3 / A_1) + (1 / A_2) + (1 / A_3) + (1 / A_4))

    ratio = Av / (D_1 * t_1)  # area ratio

    W_D = W / D_1
    t_D = t_1 / D_1
    e_D = W / (2 * D_1)

    F_tu = material_properties[material]["sigma_u"]  # ultimate
    F_ty = material_properties[material]["sigma_y"]  # yield

    F_a = F_xx  # / flanges  # forces per flange (x=a, y=b, z=c)
    F_b = F_yy / flanges
    F_c = F_zz / flanges

    M_a = F_c * l + M_xx / flanges  # moment per flange
    M_c = F_a * l + M_zz / flanges

    I_xx = (1 / 12) * t_1 * (W**3)  # moment of Inertia
    I_zz = (1 / 12) * W * (t_1**3)

    sigma_max_root = (
        F_b / (W * t_1) + (M_a * (W / 2)) / I_xx + (M_c * (t_1 / 2)) / I_zz
    )  # max stress at the root
    P_tu = (
        Kt(curve_Kt, W_D) * F_tu * t_1 * (W - D_1) * c
    )  # force it can take per failure condition
    P_Bry = K_Bry(t_D, e_D) * F_tu * D_1 * t_1
    P_ty = K_ty(curve_Kty, ratio) * F_ty * D_1 * t_1

    # if sigma_max_root < F_ty and P_tu > F_b and P_Bry > F_b and P_ty > F_c:
    #     conclusion = "pass"
    # else:
    #     conclusion = "fail"

    return F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c


def constraint_1(vars):
    W, D_1, t_1, le = vars
    F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c = load_check(
        material,
        curve_Kt,
        flanges,
        c,
        F_x,
        F_y,
        F_z,
        M_x,
        M_y,
        M_z,
        W,
        D_1,
        t_1,
        e,
        le,
        curve_Kty,
    )
    return F_ty - sigma_max_root


def constraint_2(vars):
    W, D_1, t_1, le = vars
    F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c = load_check(
        material,
        curve_Kt,
        flanges,
        c,
        F_x,
        F_y,
        F_z,
        M_x,
        M_y,
        M_z,
        W,
        D_1,
        t_1,
        e,
        le,
        curve_Kty,
    )
    return P_tu - F_b


def constraint_3(vars):
    W, D_1, t_1, le = vars
    F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c = load_check(
        material,
        curve_Kt,
        flanges,
        c,
        F_x,
        F_y,
        F_z,
        M_x,
        M_y,
        M_z,
        W,
        D_1,
        t_1,
        e,
        le,
        curve_Kty,
    )
    return P_Bry - F_b


def constraint_4(vars):
    W, D_1, t_1, le = vars
    F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c = load_check(
        material,
        curve_Kt,
        flanges,
        c,
        F_x,
        F_y,
        F_z,
        M_x,
        M_y,
        M_z,
        W,
        D_1,
        t_1,
        e,
        le,
        curve_Kty,
    )
    return P_ty - F_c
    # Forces per lug:


def constraint_5(vars):
    W, D_1, t_1, le = vars
    return W - D_1 - 0.002


def constraint_6(vars):
    W, D_1, t_1, le = vars
    return W * le * t_1 + 0.5 * np.pi * (W**2 / 4) * t_1 - np.pi * D_1**2 / 4 * t_1


def constraint_thickness_gt_width(vars):
    W, D_1, t_1, le = vars
    return W - 3 * t_1 - 0.001


def constraint_A1_positive(vars):
    W, D_1, t_1, le = vars
    return t_1 * (W / 2 + D_1 / 2 * np.sin(np.pi / 4))


def constraint_A2_positive(vars):
    W, D_1, t_1, le = vars
    return t_1 * (W / 2 - D_1 / 2)


def constraint_le_gt_D1(vars):
    W, D_1, t_1, le = vars
    return le - D_1 - 0.005


# F_x = 15000  # N
# F_y = 15000  # N
# F_z = 15000  # N
# M_x = 15  # Nm
# M_y = 0  # Nm
# M_z = 15  # Nm

# F_x = 1.597  # N
# F_y = 0.09743  # N
# F_z = 0.05389  # N
# M_x = 0.001  # Nm
# M_y = 0  # Nm
# M_z = 0  # Nm

F_x = 1594.125  # N
F_y = 1594.125  # N
F_z = 1594.125  # N
M_x = 0  # Nm
M_y = 0  # Nm
M_z = 0  # Nm
# Geometry m

W = 0.03
D_1 = 0.02
t_1 = 0.2
e = W / 2  # distance center of D_1 to the end of the flange (currently a circular end)
l = 0.1  # flange length

flanges = 2  # number of flanges


# Material Properties

material = "2014-T6"  # choose a material from the dictionary material_properties
curve_Kt = (
    "Curve 1"  # Change only the number (based on which material and kinda on geometery)
)
curve_Kty = "Curve 3"  # Change only the number


c = 0.6  # reduction due to ultimate instead or yield (tension)


a = load_check(
    material,
    curve_Kt,
    flanges,
    c,
    F_x,
    F_y,
    F_z,
    M_x,
    M_y,
    M_z,
    W,
    D_1,
    t_1,
    e,
    l,
    curve_Kty,
)


def volume(vars):
    W, D_1, t_1, le = vars
    return W * le * t_1 + 0.5 * np.pi * (W**2 / 4) * t_1 - np.pi * D_1**2 / 4 * t_1


constraints = [
    {"type": "ineq", "fun": constraint_1},
    {"type": "ineq", "fun": constraint_2},
    {"type": "ineq", "fun": constraint_3},
    {"type": "ineq", "fun": constraint_4},
    {"type": "ineq", "fun": constraint_5},
    {"type": "ineq", "fun": constraint_6},
    {"type": "ineq", "fun": constraint_thickness_gt_width},
    {"type": "ineq", "fun": constraint_A1_positive},
    {"type": "ineq", "fun": constraint_A2_positive},
    {"type": "ineq", "fun": constraint_le_gt_D1},
]
bounds = [(0.005, None), (0.005, None), (0.00001, None), (0.005, None)]
x0 = [0.0001, 0.00005, 0.0001, 0.0001]
result = minimize(volume, x0, method="SLSQP", bounds=bounds, constraints=constraints)
print("Optimization success:", result.success)
print("Minimum volume =", result.fun)
print("W = ", result.x[0])
print("D_1 = ", result.x[1])
print("t_1 = ", result.x[2])
print("le = ", result.x[3])
print(result.fun * flanges * material_properties[material]["density"])
print("F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c ")
print(
    load_check(
        material,
        curve_Kt,
        flanges,
        c,
        F_x,
        F_y,
        F_z,
        M_x,
        M_y,
        M_z,
        result.x[0],
        result.x[1],
        result.x[2],
        result.x[0] / 2,
        result.x[3],
        curve_Kty,
    )
)
print("A1", result.x[2] * (result.x[0] / 2 - result.x[1] / 2 * np.sin(np.pi / 4)))
print("A2", result.x[2] * (result.x[0] / 2 - result.x[1] / 2))
# print(a)

A_1 = result.x[2] * (result.x[0] / 2 - result.x[1] / 2 * np.sin(np.pi / 4))
A_2 = result.x[2] * (result.x[0] / 2 - result.x[1] / 2)
A_3 = A_2
A_4 = A_1
Av = 6 / (3 / A_1 + 1 / A_2 + 1 / A_3 + 1 / A_4)

ratio = Av / (result.x[1] * result.x[2])  # area ratio
print("Av = ", Av)
print("ratio = ", ratio)
print(
    f"W_D = {result.x[0] / result.x[1]}\nt_D = {result.x[3] / result.x[1]}\ne_D = {result.x[0] / (2 * result.x[1])}\nK_t = {Kt(curve_Kt, result.x[0] / result.x[1])}\nK_Bry = {K_Bry(result.x[3] / result.x[1], result.x[0] / (2 * result.x[1]))}\nK_ty = {K_ty(curve_Kty, ratio)}"
)
