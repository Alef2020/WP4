import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def load_csv_decimal_comma(filename):
    # read file as text
    with open(filename, "r") as f:
        lines = f.readlines()
    
    # replace decimal commas with decimal points
    cleaned = [line.replace(",", ".") for line in lines]
    
    # also replace semicolon delimiters with spaces
    cleaned = [line.replace(";", " ") for line in cleaned]

    # convert to numpy
    data = np.genfromtxt(cleaned, skip_header=1)
    return data

data = load_csv_decimal_comma("Kt Curves.csv")

x = data[:, 6]  # W/D
y = data[:, 7]  # Kt

degree = 4
coeffs = np.polyfit(x, y, degree)
p = np.poly1d(coeffs)

f_interp = interp1d(x, y, kind='cubic', fill_value='extrapolate')


x_smooth = np.linspace(1, 5, 200)
y_smooth = f_interp (x_smooth)

plt.scatter(x, y, label="digitized data")
plt.plot(x_smooth, y_smooth, label="polynomial fit")
plt.legend()
plt.show()
