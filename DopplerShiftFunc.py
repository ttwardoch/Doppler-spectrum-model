import numpy as np
from settings import *
from HelperFunc import *

def doppler_shift(satellite_angle_perp, satellite_angle_tang):
    """
    :param satellite_angle_perp: inclination of the satellite in degrees
    :param satellite_angle_tang: angular distance of the satellite from the highest point above the horizon (- is before achieving the highest point) in degrees
    :return: value of Doppler shift in Hz
    """
    position_of_satellite = np.array([0, SATELLITE_RADIUS*np.cos(satellite_angle_tang*np.pi/180), SATELLITE_RADIUS*np.sin(satellite_angle_tang*np.pi/180)])
    velocity_of_satellite = np.array([0, ORBITAL_SPEED*np.sin(-satellite_angle_tang*np.pi/180), ORBITAL_SPEED*np.cos(satellite_angle_tang*np.pi/180)])
    position_of_UE = np.array([EARTH_RADIUS * np.sin(satellite_angle_perp*np.pi/180), EARTH_RADIUS*np.cos(satellite_angle_perp*np.pi/180), 0])
    velocity_of_UE = np.array([EARTH_SPEED*np.cos(satellite_angle_perp*np.pi/180), EARTH_SPEED * np.sin(-satellite_angle_perp*np.pi/180),0]) * np.cos(satellite_angle_tang*np.pi/180)
    position_of_UE_minus_satellite = position_of_UE - position_of_satellite
    velocity_of_satellite_minus_UE = velocity_of_satellite - velocity_of_UE

    # Check if above horizon
    if np.dot((position_of_satellite - position_of_UE), position_of_UE) < 0:
        raise ValueError

    doppler = FREQUENCY / c * np.dot(velocity_of_satellite_minus_UE, unit_vector(position_of_UE_minus_satellite))
    return doppler

def signal_delay(satellite_angle_perp, satellite_angle_tang):
    """
    :param satellite_angle_perp: inclination of the satellite in degrees
    :param satellite_angle_tang: angular distance of the satellite from the highest point above the horizon (- is before achieving the highest point) in degrees
    :return: value of delay in ms
    """
    position_of_satellite = np.array([0, SATELLITE_RADIUS*np.cos(satellite_angle_tang*np.pi/180), SATELLITE_RADIUS*np.sin(satellite_angle_tang*np.pi/180)])
    position_of_UE = np.array([EARTH_RADIUS * np.sin(satellite_angle_perp*np.pi/180), EARTH_RADIUS*np.cos(satellite_angle_perp*np.pi/180), 0])
    position_of_UE_minus_satellite = position_of_UE - position_of_satellite

    # Check if above horizon
    if np.dot((position_of_satellite - position_of_UE), position_of_UE) < 0:
        raise ValueError

    delay = np.linalg.norm(position_of_UE_minus_satellite) / c * 1000
    #print(np.linalg.norm(position_of_UE_minus_satellite))
    return delay


def total_signal_delay(satellite_angle_perp, satellite_angle_tang):
    """
    :param satellite_angle_perp: inclination of the satellite in degrees
    :param satellite_angle_tang: angular distance of the satellite from the highest point above the horizon (- is before achieving the highest point) in degrees
    :return: value of delay in ms
    """
    return np.sqrt((satellite_angle_perp - SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE)**2 + (satellite_angle_tang - SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE)**2) * np.pi/180 * EARTH_RADIUS / c * 1000 + signal_delay(satellite_angle_perp, satellite_angle_tang)

if __name__ == '__main__':
    print(f"Orbital speed: {ORBITAL_SPEED / 1000} km/s")
    assert(doppler_shift(10, 10) == doppler_shift(-10, 10))
    assert (doppler_shift(10, -10) == -doppler_shift(10, 10))
    assert(doppler_shift(0, 0) == 0)
    assert(doppler_shift(10, 0) == 0)
    try:
        doppler_shift(90, 0)
        assert(True)
    except ValueError:
        pass

