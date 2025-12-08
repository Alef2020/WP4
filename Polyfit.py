import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def interpolation(filename, i=0, j=1):
    x_vals = []
    y_vals = []

    with open(filename, "r") as f:
        # read all lines except header
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

    # create interpolation function
    f = interp1d(x, y, kind="cubic", bounds_error = False, fill_value="extrapolate")

    return f, x, y


# Load interpolation curve + original data
f_interp, x_raw, y_raw = interpolation("K_ty Curves.csv", 4, 5)

# Smooth domain
x_smooth = np.linspace(min(x_raw), max(x_raw), 400)

# Interpolation result
y_interp = f_interp(x_smooth)

# ---- 6th-degree polynomial fit ----
coeffs_6 = np.polyfit(x_raw, y_raw, 9)   # fit polynomial
poly6 = np.polyval(coeffs_6, x_smooth)   # evaluate polynomial

# ---- Plot results ----
plt.plot(x_smooth, y_interp, label="Quadratic Interpolation", linewidth=2)
plt.plot(x_smooth, poly6, label="6th-Degree Polyfit", linestyle='--')

plt.ylim(0, 2)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Interpolation + 6th-Degree Polynomial Fit")
plt.grid(True)

plt.show()
