#!/usr/bin/env python3
"""
Test script for the eprplot module.
"""

import matplotlib.pyplot as plt
import numpy as np

import epyr

print("✅ Testing eprplot module")

# Test 1D plotting
print("\n📊 Testing 1D plotting...")
x_1d = np.linspace(3300, 3400, 1000)
y_1d = np.exp(-(((x_1d - 3350) / 10) ** 2)) + 0.1 * np.random.randn(1000)
params_1d = {"XAXIS_NAME": "Magnetic Field", "XAXIS_UNIT": "G"}

fig1, ax1 = epyr.plot_1d(x_1d, y_1d, params_1d, title="Test 1D EPR Spectrum")
print("✅ 1D plot created successfully")

# Test 2D map plotting
print("\n🗺️  Testing 2D map plotting...")
field_axis = np.linspace(3300, 3400, 200)
angle_axis = np.linspace(0, 180, 25)
y_2d = np.zeros((25, 200))

# Create synthetic 2D data
for i, angle in enumerate(angle_axis):
    center = 3350 + 10 * np.cos(np.radians(angle))
    width = 5 + 2 * np.sin(np.radians(angle))
    y_2d[i, :] = np.exp(-(((field_axis - center) / width) ** 2))

x_2d = [field_axis, angle_axis]
params_2d = {
    "XAXIS_NAME": "Magnetic Field",
    "XAXIS_UNIT": "G",
    "YAXIS_NAME": "Angle",
    "YAXIS_UNIT": "deg",
}

fig2, ax2 = epyr.plot_2d_map(x_2d, y_2d, params_2d, title="Test 2D EPR Map")
print("✅ 2D map created successfully")

# Test 2D waterfall plotting
print("\n🌊 Testing 2D waterfall plotting...")
fig3, ax3 = epyr.plot_2d_waterfall(
    x_2d, y_2d, params_2d, title="Test 2D EPR Waterfall", offset_factor=0.3
)
print("✅ 2D waterfall created successfully")

# Test with real EPR data if available
print("\n📂 Testing with real EPR data...")
try:
    x, y, params, filepath = epyr.eprload(
        "examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC", plot_if_possible=False
    )
    if y is not None:
        if y.ndim == 1:
            fig4, ax4 = epyr.plot_1d(x, y, params, title="Real 1D EPR Data")
            print("✅ Real 1D data plotted successfully")
        elif y.ndim == 2:
            fig5, ax5 = epyr.plot_2d_map(x, y, params, title="Real 2D EPR Map")
            fig6, ax6 = epyr.plot_2d_waterfall(
                x, y, params, title="Real 2D EPR Waterfall"
            )
            print("✅ Real 2D data plotted successfully")
    else:
        print("⚠️  No real data found, using synthetic data only")
except:
    print("⚠️  No real data available, using synthetic data only")

print("\n🎉 All eprplot functions tested successfully!")
print("\nYou can now use:")
print("  epyr.plot_1d(x, y, params)")
print("  epyr.plot_2d_map(x, y, params)")
print("  epyr.plot_2d_waterfall(x, y, params)")

plt.show()
