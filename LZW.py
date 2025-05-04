'''import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D

# Array Factor for ULA in 3D
def array_factor_3d(N, d, theta, phi):
    k = 2 * np.pi  # assume λ = 1
    psi = k * d * np.cos(theta)
    with np.errstate(divide='ignore', invalid='ignore'):
        AF = np.sin(N * psi / 2) / (N * np.sin(psi / 2))
        AF = np.nan_to_num(AF)
    AF = np.abs(AF)
    AF = AF / np.max(AF)
    return AF

def plot_3d_pattern():
    try:
        N = int(num_elements.get())
        d = float(element_spacing.get())
    except ValueError:
        return

    # Spherical coordinate grid
    theta = np.linspace(0, np.pi, 180)   # 0 to 180°
    phi = np.linspace(0, 2 * np.pi, 360) # 0 to 360°
    theta, phi = np.meshgrid(theta, phi)

    # Radiation intensity (only depends on theta for ULA along z-axis)
    AF = array_factor_3d(N, d, theta, phi)
    r = AF

    # Convert to Cartesian coordinates for 3D plotting
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    ax.clear()
    ax.plot_surface(x, y, z, cmap='viridis', edgecolor='k', linewidth=0.1, alpha=0.9)
    ax.set_title("3D Radiation Pattern (ULA)", pad=20)
    ax.set_box_aspect([1,1,1])
    canvas.draw()

# GUI Setup
root = tk.Tk()
root.title("3D Antenna Array Radiation Pattern")

frame = ttk.Frame(root, padding="10")
frame.grid(row=0, column=0, sticky="n")

ttk.Label(frame, text="Number of Elements (N):").grid(row=0, column=0, sticky='w')
num_elements = ttk.Entry(frame)
num_elements.insert(0, "8")
num_elements.grid(row=0, column=1)

ttk.Label(frame, text="Element Spacing (d in λ):").grid(row=1, column=0, sticky='w')
element_spacing = ttk.Entry(frame)
element_spacing.insert(0, "0.5")
element_spacing.grid(row=1, column=1)

plot_button = ttk.Button(frame, text="Plot 3D Pattern", command=plot_3d_pattern)
plot_button.grid(row=2, column=0, columnspan=2, pady=10)

# 3D Plot
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection='3d')
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=0, column=1, rowspan=3)

root.mainloop()'''

import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import messagebox

def compute_pattern():
    try:
        n = int(entry_n.get())
        d = float(entry_d.get())
        beta = float(entry_beta.get()) * np.pi / 180  # convert to radians
        amps = [float(a) for a in entry_amps.get().split(',')]

        if len(amps) != n:
            raise ValueError("Number of amplitudes must match number of elements.")

        theta = np.linspace(0, 2 * np.pi, 1000)
        k = 2 * np.pi  # wavelength λ assumed to be 1, so k = 2π/λ = 2π

        AF = np.zeros_like(theta, dtype=complex)
        for i in range(n):
            AF += amps[i] * np.exp(1j * (i * k * d * np.cos(theta) + i * beta))

        AF_normalized = np.abs(AF) / np.max(np.abs(AF))

        plt.figure(figsize=(6, 6))
        ax = plt.subplot(111, polar=True)
        ax.plot(theta, AF_normalized)
        ax.set_title("Radiation Pattern (Normalized)", va='bottom')
        plt.show()

    except Exception as e:
        messagebox.showerror("Input Error", str(e))

# GUI Setup
root = Tk()
root.title("Linear Antenna Array Radiation Pattern Simulator")

Label(root, text="Number of Elements (n):").grid(row=0, column=0, sticky=W)
entry_n = Entry(root)
entry_n.grid(row=0, column=1)

Label(root, text="Inter-element Distance (d in λ):").grid(row=1, column=0, sticky=W)
entry_d = Entry(root)
entry_d.grid(row=1, column=1)

Label(root, text="Progressive Phase Shift (β in degrees):").grid(row=2, column=0, sticky=W)
entry_beta = Entry(root)
entry_beta.grid(row=2, column=1)

Label(root, text="Amplitudes (comma-separated):").grid(row=3, column=0, sticky=W)
entry_amps = Entry(root)
entry_amps.grid(row=3, column=1)

Button(root, text="Simulate", command=compute_pattern).grid(row=4, column=0, columnspan=2, pady=10)

root.mainloop()
Hello
Hi
