import numpy as np

from constants import *


# Set the carrier frequency
FREQUENCY = 30 # GHz

# Doppler shift within a beam
PATCH_DIAMETER = 100000 # m
SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE = 0  # degrees
SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE = 0  # degrees
SATELLITE_HEIGHT = 600000 # in m
n = 100  # resolution of the coloured patch n x n

# Delays
X_SHIFT = -0.0  # How far from the middle is the UE in ANGULAR_RADIUS [-1, 1]
Y_SHIFT = 0.0  # How far from the middle is the UE in ANGULAR_RADIUS [-1, 1]

# Scatterers simulation
"""
constant_for_plots = "RX_width"  # "RX_width" or "RX_angle" determines a type of folder created and additional data
RX_widths = [5, 10, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360]  # degrees
RX_angles = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5, 360]
"""
#"""
constant_for_plots = "RX_angle"  # "RX_width" or "RX_angle" determines a type of folder created and additional data
RX_widths = [1, 2, 3, 5, 10, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360]  # degrees
RX_angles = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5, 360]
#"""
n_bins = 101  # granularity of bar diagrams
iterations = 1000  # number of random points

RX_widths = [1, 2, 3, 5, 10, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360]
RX_angles = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5, 360]

RX_widths_export = [360, ]  # degrees
RX_angles_export = [90, ]

# Type of randomness
random_type = "gaussian"  # "uniform" or "radial" or "gaussian"
STD_DEVIATION = 0.2 #  as a fraction of the radius, only for Gaussian

# Export of data
EXPORT = True
RAYLEIGH = True
if EXPORT:
    RAYLEIGH = EXPORT
    RX_widths_export = [360, ]  # degrees
    RX_angles_export = [90, ]

    RX_widths = RX_widths_export
    RX_angles = RX_angles_export

#################################################################################################
# Helper constants for the beam
PATCH_RADIUS = PATCH_DIAMETER / 2
ANGULAR_WIDTH_OF_PATCH = PATCH_DIAMETER/EARTH_RADIUS * 180/np.pi
X_SHIFT *= ANGULAR_WIDTH_OF_PATCH/2
Y_SHIFT *= ANGULAR_WIDTH_OF_PATCH/2

# Helper constants for Doppler shift estimation
SATELLITE_RADIUS = EARTH_RADIUS + SATELLITE_HEIGHT
ORBITAL_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS) # m/s
ANGULAR_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS) / SATELLITE_RADIUS
EARTH_SPEED = 2*np.pi*EARTH_RADIUS / 24 / 3600  # speed of the equator in m/s

# Conversion of frequency GHz -> Hz
FREQUENCY *= 1000000000