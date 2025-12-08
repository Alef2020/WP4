import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import time
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------
# material_properties (your original)
# ---------------------------
material_properties = {
    "2014-T6": {
        "density": 2800,
        "E": 72.4 * (10**9),
        "sigma_y": 414 * (10**6),
        "sigma_u": 483 * (10**6),
        "shear": 290 * (10**6),
    },
    "7075-T6": {
        "density": 2810,
        "E": 71.7 * (10**9),
        "sigma_y": 503 * (10**6),
        "sigma_u": 572 * (10**6),
        "shear": 331 * (10**6),
    },
    "2024-T2": {
        "density": 2780,
        "E": 73.1 * (10**9),
        "sigma_y": 324 * (10**6),
        "sigma_u": 469 * (10**6),
        "shear": 283 * (10**6),
    },
    "2024-T3": {
        "density": 2780,
        "E": 73.1 * (10**9),
        "sigma_y": 345 * (10**6),
        "sigma_u": 483 * (10**6),
        "shear": 283 * (10**6),
    },
    "2024-T4": {
        "density": 2780,
        "E": 71.0 * (10**9),
        "sigma_y": 324 * (10**6),
        "sigma_u": 469 * (10**6),
        "shear": 290 * (10**6),
    },
    "AZ1916-T6": {
        "density": 1810,
        "E": 44.8 * (10**9),
        "sigma_y": 145 * (10**6),
        "sigma_u": 275 * (10**6),
        "shear": 110 * (10**6),
    },
    "356-T6": {
        "density": 2680,
        "E": 72.4 * (10**9),
        "sigma_y": 138 * (10**6),
        "sigma_u": 207 * (10**6),
        "shear": 180 * (10**6),
    },
    "4130-steel": {
        "density": 7850,
        "E": 205 * (10**9),
        "sigma_y": 435 * (10**6),
        "sigma_u": 670 * (10**6),
        "shear": 338 * (10**6),
    },
    "8630-steel": {
        "density": 7850,
        "E": 200 * (10**9),
        "sigma_y": 550 * (10**6),
        "sigma_u": 620 * (10**6),
        "shear": 338 * (10**6),
    },
}

# ---------------------------
# helper functions (polyfit interpolation as before)
# ---------------------------
def interpolation(filename, i=0, j=1):
    x_vals = []
    y_vals = []
    with open(filename, "r") as f:
        lines = f.readlines()[1:]
    for line in lines:
        line = line.replace(",", ".").replace(";", " ")
        parts = line.strip().split()
        if len(parts) <= max(i, j):
            continue
        try:
            x = float(parts[i])
            y = float(parts[j])
            x_vals.append(x)
            y_vals.append(y)
        except ValueError:
            continue
    x = np.array(x_vals)
    y = np.array(y_vals)
    coeffs = np.polyfit(x, y, min(9, max(1, len(x) - 1)))
    return coeffs


def Kt(curve, W_D):
    if curve == "Curve 1":
        f_interp = interpolation("Kt Curves.csv", 0, 1)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 2":
        f_interp = interpolation("Kt Curves.csv", 2, 3)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 3":
        f_interp = interpolation("Kt Curves.csv", 4, 5)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 4":
        f_interp = interpolation("Kt Curves.csv", 6, 7)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 5":
        f_interp = interpolation("Kt Curves.csv", 8, 9)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 6":
        f_interp = interpolation("Kt Curves.csv", 10, 11)
        K_t = np.polyval(f_interp, W_D)
    elif curve == "Curve 7":
        f_interp = interpolation("Kt Curves.csv", 12, 13)
        K_t = np.polyval(f_interp, W_D)
    else:
        raise ValueError("Unknown Kt curve")
    return K_t


def K_Bry(t_D, e_D):
    if t_D <= 0.07:
        f_interp = interpolation("K_bry Curves.csv", 0, 1)
    elif t_D <= 0.09:
        f_interp = interpolation("K_bry Curves.csv", 2, 3)
    elif t_D <= 0.11:
        f_interp = interpolation("K_bry Curves.csv", 4, 5)
    elif t_D <= 0.135:
        f_interp = interpolation("K_bry Curves.csv", 6, 7)
    elif t_D <= 0.175:
        f_interp = interpolation("K_bry Curves.csv", 8, 9)
    elif t_D <= 0.25:
        f_interp = interpolation("K_bry Curves.csv", 10, 11)
    elif t_D <= 0.35:
        f_interp = interpolation("K_bry Curves.csv", 12, 13)
    elif t_D <= 0.5:
        f_interp = interpolation("K_bry Curves.csv", 14, 15)
    else:
        f_interp = interpolation("K_bry Curves.csv", 16, 17)
    return np.polyval(f_interp, e_D)


def K_ty(curve, Av_Abr):
    mapping = {
        "Curve 1": (0, 1),
        "Curve 2": (2, 3),
        "Curve 3": (4, 5),
        "Curve 4": (6, 7),
        "Curve 5": (8, 9),
        "Curve 6": (10, 11),
        "Curve 7": (12, 13),
        "Curve 8": (14, 15),
        "Curve 9": (16, 17),
        "Curve 10": (18, 19),
        "Curve 11": (20, 21),
        "Curve 12": (22, 23),
        "Curve 13": (24, 25),
        "Curve 14": (26, 27),
    }
    if curve not in mapping:
        raise ValueError("Unknown K_ty curve")
    i, j = mapping[curve]
    coeffs = interpolation("K_ty Curves.csv", i, j)
    return np.polyval(coeffs, Av_Abr)


# ---------------------------
# load_check (unchanged logic)
# ---------------------------
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
    A_1 = t_1 * (W / 2 - (D_1 / 2) * np.sin(np.pi / 4))
    A_2 = t_1 * (W / 2 - D_1 / 2)
    A_3 = A_2
    A_4 = A_1
    if min(A_1, A_2, A_3, A_4) <= 0:
        Av = 0.0
    else:
        Av = 6 / ((3 / A_1) + (1 / A_2) + (1 / A_3) + (1 / A_4))
    ratio = Av / (D_1 * t_1) if D_1 * t_1 != 0 else 0

    W_D = W / D_1 if D_1 != 0 else 1e6
    t_D = t_1 / D_1 if D_1 != 0 else 1e6
    e_D = W / (2 * D_1) if D_1 != 0 else 1e6

    F_tu = material_properties[material]["sigma_u"]
    F_ty = material_properties[material]["sigma_y"]

    F_a = F_xx
    F_b = F_yy / flanges
    F_c = F_zz / flanges

    M_a = F_c * l + M_xx / flanges
    M_c = F_a * l + M_zz / flanges

    I_xx = (1 / 12) * t_1 * (W**3)
    I_zz = (1 / 12) * W * (t_1**3)

    sigma_max_root = 1e12
    try:
        sigma_max_root = (
            F_b / (W * t_1) + (M_a * (W / 2)) / I_xx + (M_c * (t_1 / 2)) / I_zz
        )
    except Exception:
        sigma_max_root = 1e12

    try:
        P_tu = Kt(curve_Kt, W_D) * F_tu * t_1 * (W - D_1) * c
    except Exception:
        P_tu = -1e12
    try:
        P_Bry = K_Bry(t_D, e_D) * F_ty * D_1 * t_1
    except Exception:
        P_Bry = -1e12
    try:
        P_ty = K_ty(curve_Kty, ratio) * F_ty * D_1 * t_1
    except Exception:
        P_ty = -1e12

    F_shear_max = (F_b**2 + F_c**2) ** (1 / 2)
    sigma_bolt_shear = (
        F_shear_max / (2 * np.pi * ((D_1 / 2) ** 2)) if D_1 > 0 else 1e12
    )

    return F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear


# ---------------------------
# problem inputs
# ---------------------------
F_x = 490.5  # N
F_y = 490.5  # N
F_z = 1373.4  # N
M_x = 0  # Nm
M_y = 0  # Nm
M_z = 0  # Nm

flanges = 2
material = "2014-T6"
curve_Kt = "Curve 1"
curve_Kty = "Curve 3"
c = 0.6

# bounds
W_bounds = (0.01, 0.05)     # [m] sweep range for W
D1_bounds = (0.005, 0.025)  # [m] sweep range for D_1
t1_bounds = (0.001, 0.05)   # [m] thickness
le_bounds = (0.01, 0.2)     # [m] flange length

# grid resolution
n_W = 40
n_D1 = 40
n_t1 = 40

# ---------------------------
# volume function
# ---------------------------
def volume_from_params(W, D1, t1, le):
    return W * le * t1 + 0.5 * np.pi * (W**2 / 4) * t1 - np.pi * (D1**2 / 4) * t1


# ===========================
# SWEEP 1: Surface vs (W, D1)
# optimize over (t1, le) for each (W,D1)
# ===========================
W_grid = np.linspace(W_bounds[0], W_bounds[1], n_W)
D1_grid = np.linspace(D1_bounds[0], D1_bounds[1], n_D1)

# Use indexing='ij' so shapes are (n_W, n_D1) and match loop ordering
W_mesh, D1_mesh = np.meshgrid(W_grid, D1_grid, indexing='ij')

mass_grid_WD = np.full((n_W, n_D1), np.nan)
feasible_mask_WD = np.zeros((n_W, n_D1), dtype=bool)
opt_t1_WD = np.full((n_W, n_D1), np.nan)
opt_le_WD = np.full((n_W, n_D1), np.nan)

start_time = time.time()
total = n_W * n_D1
count = 0

for i in range(n_W):
    for j in range(n_D1):
        count += 1
        W_val = W_grid[i]
        D1_val = D1_grid[j]

        # quick infeasible checks
        if D1_val >= W_val - 1e-8:
            continue

        def inner_obj(x):
            t1, le = x
            return volume_from_params(W_val, D1_val, t1, le)

        # constraints for fixed W_val, D1_val
        def c1(x):
            t1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, D1_val, t1, W_val/2, le, curve_Kty
            )
            return F_ty - sigma_max_root

        def c2(x):
            t1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, D1_val, t1, W_val/2, le, curve_Kty
            )
            return P_tu - F_b

        def c3(x):
            t1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, D1_val, t1, W_val/2, le, curve_Kty
            )
            return P_Bry - F_b

        def c4(x):
            t1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, D1_val, t1, W_val/2, le, curve_Kty
            )
            return P_ty - F_c

        def c5(x):
            t1, le = x
            return W_val - D1_val - 0.002

        def c6(x):
            t1, le = x
            return W_val * le * t1 + 0.5 * np.pi * (W_val**2 / 4) * t1 - np.pi * (D1_val**2 / 4) * t1

        def c_thickness_gt_width(x):
            t1, le = x
            return W_val - t1 - 0.001

        def c_A1_positive(x):
            t1, le = x
            return (W_val / 2 - (D1_val / 2) * np.sin(np.pi / 4)) * t1

        def c_A2_positive(x):
            t1, le = x
            return (W_val / 2 - D1_val / 2) * t1

        def c_le_gt_D1(x):
            t1, le = x
            return le - D1_val - 0.005

        def c_sigma_bolt_shear(x):
            t1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, D1_val, x[0], W_val/2, x[1], curve_Kty
            )
            return material_properties[material]["shear"] * c - sigma_bolt_shear

        cons = [
            {"type": "ineq", "fun": c1},
            {"type": "ineq", "fun": c2},
            {"type": "ineq", "fun": c3},
            {"type": "ineq", "fun": c4},
            {"type": "ineq", "fun": c5},
            {"type": "ineq", "fun": c6},
            {"type": "ineq", "fun": c_thickness_gt_width},
            {"type": "ineq", "fun": c_A1_positive},
            {"type": "ineq", "fun": c_A2_positive},
            {"type": "ineq", "fun": c_le_gt_D1},
            {"type": "ineq", "fun": c_sigma_bolt_shear},
        ]

        x0 = [max(t1_bounds[0], 0.005), max(le_bounds[0], D1_val + 0.01)]
        bnds = [t1_bounds, le_bounds]

        try:
            res = minimize(inner_obj, x0, method="SLSQP", bounds=bnds, constraints=cons, options={"maxiter": 300, "ftol": 1e-6})
        except Exception:
            res = None

        if res is not None and res.success:
            vol = inner_obj(res.x)
            mass = vol * flanges * material_properties[material]["density"]
            mass_grid_WD[i, j] = mass
            feasible_mask_WD[i, j] = True
            opt_t1_WD[i, j] = res.x[0]
            opt_le_WD[i, j] = res.x[1]
        else:
            mass_grid_WD[i, j] = np.nan

        if count % 200 == 0 or count == total:
            elapsed = time.time() - start_time
            print(f"[W,D1 sweep] Processed {count}/{total} points, elapsed {elapsed:.1f}s")

# ===========================
# SWEEP 2: Surface vs (W, t1)
# optimize over (D1, le) for each (W,t1)
# ===========================
W_grid2 = np.linspace(W_bounds[0], W_bounds[1], n_W)
t1_grid = np.linspace(t1_bounds[0], t1_bounds[1], n_t1)

W_mesh2, t1_mesh = np.meshgrid(W_grid2, t1_grid, indexing='ij')

mass_grid_Wt1 = np.full((n_W, n_t1), np.nan)
feasible_mask_Wt1 = np.zeros((n_W, n_t1), dtype=bool)
opt_D1_Wt1 = np.full((n_W, n_t1), np.nan)
opt_le_Wt1 = np.full((n_W, n_t1), np.nan)

start_time2 = time.time()
total2 = n_W * n_t1
count2 = 0

for i in range(n_W):
    for j in range(n_t1):
        count2 += 1
        W_val = W_grid2[i]
        t1_val = t1_grid[j]

        # quick infeasible checks
        if t1_val >= W_val - 1e-8:
            # thickness cannot exceed width (constraint)
            continue

        def inner_obj2(x):
            D1, le = x
            return volume_from_params(W_val, D1, t1_val, le)

        # constraints adapted for fixed W_val, t1_val
        def c1b(x):
            D1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, x[0], t1_val, W_val/2, x[1], curve_Kty
            )
            return F_ty - sigma_max_root

        def c2b(x):
            D1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, x[0], t1_val, W_val/2, x[1], curve_Kty
            )
            return P_tu - F_b

        def c3b(x):
            D1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, x[0], t1_val, W_val/2, x[1], curve_Kty
            )
            return P_Bry - F_b

        def c4b(x):
            D1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, x[0], t1_val, W_val/2, x[1], curve_Kty
            )
            return P_ty - F_c

        def c5b(x):
            D1, le = x
            return W_val - x[0] - 0.002

        def c6b(x):
            D1, le = x
            return W_val * le * t1_val + 0.5 * np.pi * (W_val**2 / 4) * t1_val - np.pi * (x[0]**2 / 4) * t1_val

        def c_thickness_gt_width_b(x):
            D1, le = x
            return W_val - t1_val - 0.001

        def c_A1_positive_b(x):
            D1, le = x
            return (W_val / 2 - (x[0] / 2) * np.sin(np.pi / 4)) * t1_val

        def c_A2_positive_b(x):
            D1, le = x
            return (W_val / 2 - x[0] / 2) * t1_val

        def c_le_gt_D1_b(x):
            D1, le = x
            return le - x[0] - 0.005

        def c_sigma_bolt_shear_b(x):
            D1, le = x
            F_ty, sigma_max_root, P_tu, P_Bry, F_b, P_ty, F_c, sigma_bolt_shear = load_check(
                material, curve_Kt, flanges, c,
                F_x, F_y, F_z, M_x, M_y, M_z,
                W_val, x[0], t1_val, W_val/2, x[1], curve_Kty
            )
            return material_properties[material]["shear"] * c - sigma_bolt_shear

        cons_b = [
            {"type": "ineq", "fun": c1b},
            {"type": "ineq", "fun": c2b},
            {"type": "ineq", "fun": c3b},
            {"type": "ineq", "fun": c4b},
            {"type": "ineq", "fun": c5b},
            {"type": "ineq", "fun": c6b},
            {"type": "ineq", "fun": c_thickness_gt_width_b},
            {"type": "ineq", "fun": c_A1_positive_b},
            {"type": "ineq", "fun": c_A2_positive_b},
            {"type": "ineq", "fun": c_le_gt_D1_b},
            {"type": "ineq", "fun": c_sigma_bolt_shear_b},
        ]

        x0b = [max(D1_bounds[0], 0.01 * W_val), max(le_bounds[0], D1_bounds[0] + 0.02)]
        bndsb = [D1_bounds, le_bounds]

        try:
            res2 = minimize(inner_obj2, x0b, method="SLSQP", bounds=bndsb, constraints=cons_b, options={"maxiter": 300, "ftol": 1e-6})
        except Exception:
            res2 = None

        if res2 is not None and res2.success:
            vol2 = inner_obj2(res2.x)
            mass2 = vol2 * flanges * material_properties[material]["density"]
            mass_grid_Wt1[i, j] = mass2
            feasible_mask_Wt1[i, j] = True
            opt_D1_Wt1[i, j] = res2.x[0]
            opt_le_Wt1[i, j] = res2.x[1]
        else:
            mass_grid_Wt1[i, j] = np.nan

        if count2 % 200 == 0 or count2 == total2:
            elapsed = time.time() - start_time2
            print(f"[W,t1 sweep] Processed {count2}/{total2} points, elapsed {elapsed:.1f}s")

# -------------------
# Plot Surface (W, D1)
# -------------------
fig1 = plt.figure(figsize=(11, 7))
ax1 = fig1.add_subplot(111, projection='3d')
W_plot = W_mesh
D1_plot = D1_mesh
Z1 = mass_grid_WD
Z1_masked = np.ma.masked_invalid(Z1)

surf1 = ax1.plot_surface(W_plot, D1_plot, Z1_masked, cmap='viridis', edgecolor='none', linewidth=0, antialiased=True)
ax1.set_xlabel('W [m]')
ax1.set_ylabel('D_1 [m]')
ax1.set_zlabel('Mass [kg]')
ax1.set_title(f"Mass surface vs (W, D1) (material={material}, flanges={flanges})")
fig1.colorbar(surf1, shrink=0.6, aspect=10, label='Mass [kg]')
ax1.contourf(W_plot, D1_plot, Z1_masked, zdir='z', offset=np.nanmin(Z1_masked)-0.1*np.nanmean(Z1_masked), cmap='viridis', alpha=0.6)
plt.tight_layout()

# -------------------
# Plot Surface (W, t1)
# -------------------
fig2 = plt.figure(figsize=(11, 7))
ax2 = fig2.add_subplot(111, projection='3d')
W_plot2 = W_mesh2
t1_plot = t1_mesh
Z2 = mass_grid_Wt1
Z2_masked = np.ma.masked_invalid(Z2)

surf2 = ax2.plot_surface(W_plot2, t1_plot, Z2_masked, cmap='plasma', edgecolor='none', linewidth=0, antialiased=True)
ax2.set_xlabel('W [m]')
ax2.set_ylabel('t_1 [m]')
ax2.set_zlabel('Mass [kg]')
ax2.set_title(f"Mass surface vs (W, t1) (material={material}, flanges={flanges})")
fig2.colorbar(surf2, shrink=0.6, aspect=10, label='Mass [kg]')
ax2.contourf(W_plot2, t1_plot, Z2_masked, zdir='z', offset=np.nanmin(Z2_masked)-0.1*np.nanmean(Z2_masked), cmap='plasma', alpha=0.6)
plt.tight_layout()

plt.show()

# -------------------
# Save results (including last point values)
# -------------------
np.savez(
    "mass_surfaces_full_data.npz",
    W_grid=W_grid,
    D1_grid=D1_grid,
    t1_grid=t1_grid,
    mass_grid_WD=mass_grid_WD,
    feasible_mask_WD=feasible_mask_WD,
    opt_t1_WD=opt_t1_WD,
    opt_le_WD=opt_le_WD,
    mass_grid_Wt1=mass_grid_Wt1,
    feasible_mask_Wt1=feasible_mask_Wt1,
    opt_D1_Wt1=opt_D1_Wt1,
    opt_le_Wt1=opt_le_Wt1,
)

print("Done.")
print("WD valid fraction:", np.count_nonzero(~np.isnan(mass_grid_WD)) / mass_grid_WD.size)
print("Wt1 valid fraction:", np.count_nonzero(~np.isnan(mass_grid_Wt1)) / mass_grid_Wt1.size)
# Example: print last stored values (explicitly show last grid point data)
last_i = n_W - 1
last_j_D1 = n_D1 - 1
last_j_t1 = n_t1 - 1
print("--- Last stored entries ---")
print("Last (W,D1) point: W=", W_grid[last_i], " D1=", D1_grid[last_j_D1])
print("mass (W,D1) at last:", mass_grid_WD[last_i, last_j_D1])
print("opt t1, le (W,D1):", opt_t1_WD[last_i, last_j_D1], opt_le_WD[last_i, last_j_D1])
print("Last (W,t1) point: W=", W_grid2[last_i], " t1=", t1_grid[last_j_t1])
print("mass (W,t1) at last:", mass_grid_Wt1[last_i, last_j_t1])
print("opt D1, le (W,t1):", opt_D1_Wt1[last_i, last_j_t1], opt_le_Wt1[last_i, last_j_t1])
