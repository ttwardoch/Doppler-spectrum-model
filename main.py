import random
import os

import matplotlib.pyplot as plt
from tqdm import tqdm

from settings import *
from DopplerShiftFunc import *


# Create a grid for a patch of the Earth within the beam, assume it is flat
X = np.linspace(0, ANGULAR_WIDTH_OF_PATCH, num=n) - ANGULAR_WIDTH_OF_PATCH / 2
Y = np.linspace(0, ANGULAR_WIDTH_OF_PATCH, num=n) - ANGULAR_WIDTH_OF_PATCH / 2


# Tables for collecting variables
xs = []
ys = []
shifts_all = []
delays_all = []
delays_all_total = []
for x in tqdm(X):
    for y in Y:
        if np.sqrt(x**2 + y**2) <= ANGULAR_WIDTH_OF_PATCH/2:
            xs.append(x)
            ys.append(y)
            shifts_all.append(doppler_shift(SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE + x, SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE + y))
            delays_all.append(signal_delay(SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE + x, SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE + y))
            delays_all_total.append(np.sqrt((x-X_SHIFT)**2 + (y-Y_SHIFT)**2) * np.pi/180 * EARTH_RADIUS / c * 1000 + signal_delay(SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE + x, SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE + y))
        else:
            xs.append(x)
            ys.append(y)
            shifts_all.append(-1000000000)
            delays_all.append(-1000000000)
            delays_all_total.append(-1000000000)

# Get rid of shifts outside of the patch for min and max
shifts_valid = [shift for shift in shifts_all if shift > -1000000000]
delays_valid = [dela for dela in delays_all if dela > -1000000000]
delays_total_valid = [dela for dela in delays_all_total if dela > -1000000000]
print(f"Doppler spread: {max(shifts_valid) - min(shifts_valid)}")


if random_type == "uniform":
    main_folder = "Plots_uniform"
elif random_type == "radial":
    main_folder = "Plots_radial"

elif random_type == "gaussian":
    main_folder = f"Plots_gaussian_0p{int(STD_DEVIATION*10)}"
else:
    raise ValueError

if not os.path.exists(f"{main_folder}"):
    os.makedirs(f"{main_folder}")

# Create a folder for the plots
folder = f"SAN_{SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE}_{SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE}_PATCH_{PATCH_DIAMETER//1000}"
if not os.path.exists(f"{main_folder}/{folder}"):
    os.makedirs(f"{main_folder}/{folder}")

# Plot the Doppler shift distribution in the beam
plot_2D(shifts_all, xs, ys, n, crange=[min(shifts_valid), max(shifts_valid)], title=f"Doppler shifts in a beam {PATCH_DIAMETER}km wide\n orbit inclination:{SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE}, angle to highest elevation: {SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE}\n max Doppler spread: {int(max(shifts_valid) - min(shifts_valid))} Hz\n frequency: {int(FREQUENCY/1e9)}GHz")
plt.savefig(f"{main_folder}/{folder}/Doppler.png")

plot_2D(shifts_all, xs, ys, n, crange=[min(shifts_valid), max(shifts_valid)], title=f"Doppler shift distribution in a beam [Hz]")
plt.xlabel("x-axis [km]")
plt.ylabel("y-axis [km]")
plt.tight_layout()
plt.savefig(f"{main_folder}/{folder}/Doppler_ready.png")

# plot delays of points
plot_2D(delays_all, xs, ys, n, crange=[min(delays_valid), max(delays_valid)], title=f"Delays [ms] in a beam {PATCH_DIAMETER}km wide\n orbit inclination:{SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE}, angle to highest elevation: {SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE}\n max Delay spread: {max(delays_valid) - min(delays_valid)} Hz\n frequency: {int(FREQUENCY/1e9)}GHz")
plt.savefig(f"{main_folder}/{folder}/Delays.png")

plot_2D(delays_all, xs, ys, n, crange=[min(delays_valid), max(delays_valid)], title=f"Delays in a beam [ms]")
plt.xlabel("x-axis [km]")
plt.ylabel("y-axis [km]")
plt.tight_layout()
plt.savefig(f"{main_folder}/{folder}/Delays_ready.png")

# plot delays of points with respect to UE
plot_2D(delays_all_total, xs, ys, n, crange=[min(delays_total_valid), max(delays_total_valid)], title=f"Delays [ms] in a beam {PATCH_DIAMETER}km wide\n orbit inclination:{SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE}, angle to highest elevation: {SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE}\n max Delay spread: {max(delays_total_valid) - min(delays_total_valid)} Hz\n frequency: {int(FREQUENCY/1e9)}GHz")
plt.savefig(f"{main_folder}/{folder}/Delays_for_UE.png")

plot_2D(delays_all_total, xs, ys, n, crange=[min(delays_total_valid), max(delays_total_valid)], title=f"Delays in a beam for UE [ms]")
plt.xlabel("x-axis [km]")
plt.ylabel("y-axis [km]")
plt.tight_layout()
plt.savefig(f"{main_folder}/{folder}/Delays_for_UE_ready.png")
#plt.show()

#raise ValueError


if constant_for_plots == "RX_width":
    for RX_width in RX_widths:

        subfolder = f"RXwidth_{RX_width}"
        if not os.path.exists(f"{main_folder}/{folder}/{subfolder}"):
            os.makedirs(f"{main_folder}/{folder}/{subfolder}")

        for RX_angle in tqdm(RX_angles):
            # Jakes model but for scaterrers
            MAX_DOPPLER = max(shifts_valid)
            MIN_DOPPLER = min(shifts_valid)
            DELTA = (MAX_DOPPLER - MIN_DOPPLER)/n_bins

            bins = np.zeros(n_bins, dtype=np.complex128)
            bins_labels = np.linspace(MIN_DOPPLER, MAX_DOPPLER, num=n_bins+1) + DELTA/2
            bins_labels = bins_labels[:-1].copy()

            for _ in range(iterations):
                angle_perp, angle_tang = random_angles(random_type, ANGULAR_WIDTH_OF_PATCH)
                # Check if within circle
                if np.sqrt(angle_tang**2 + angle_perp**2) < ANGULAR_WIDTH_OF_PATCH/2:
                    #  Limit receiver angle
                    if cos_of_angle_between(unit_vector_at_angle(RX_angle), [angle_tang, angle_perp]) > receiver_angle(RX_width):
                        angle_perp += SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE
                        angle_tang += SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE
                        doppler = doppler_shift(angle_perp, angle_tang)
                        bin_target = int((doppler-MIN_DOPPLER) // DELTA)
                        if bin_target >= n_bins:
                            bin_target = n_bins-1
                        if bin_target <= -1:
                            bin_target = 0
                        if RAYLEIGH:
                            bins[bin_target] += random_unit_circle_complex_value()/total_signal_delay(angle_perp, angle_tang)/total_signal_delay(angle_perp, angle_tang)*4
                        else:
                            bins[bin_target] += 1/total_signal_delay(angle_perp, angle_tang)/total_signal_delay(angle_perp, angle_tang)*4

            bins = np.abs(bins)
            if not RAYLEIGH:
                bins = np.sqrt(bins)

            # Plot the combined figure
            plt.figure(0)
            plt.plot(bins_labels, bins)

            # Plot the bar plots and save to folder
            plt.figure(random.randint(1, 10000000000))
            plt.bar(bins_labels, bins, width=DELTA)
            plt.xlabel("Frequency[Hz]")
            plt.ylabel("Amplitude")
            plt.title(f"Doppler spectrum for RXwidth={RX_width} and RXangle={int(RX_angle)}")
            plt.title(f"Doppler spectrum")

            file = f"RXwidth_{RX_width}_Angle_{int(RX_angle)}.png"
            plt.savefig(f"{main_folder}/{folder}/{subfolder}/{file}")
            #plt.show(block=False)
            plt.close()

        plt.figure(0)
        plt.xlabel("Frequency[Hz]")
        plt.ylabel("Amplitude")
        plt.title(f"Doppler spectra for RXwidth = {RX_width}")
        plt.savefig(f"{main_folder}/{folder}/{subfolder}/RXwidth_{RX_width}.png")
        plt.clf()
        print(f"Successfully saved to: {main_folder}/" + folder + "/" + subfolder)
        #plt.show()

if constant_for_plots == "RX_angle":
    for RX_angle in RX_angles:
        print(f"Working on RX_angle = {RX_angle}")
        subfolder = f"RXangle_{int(RX_angle)}"
        if not os.path.exists(f"{main_folder}/{folder}/{subfolder}"):
            os.makedirs(f"{main_folder}/{folder}/{subfolder}")

        spreads = []
        for RX_width in tqdm(RX_widths):
            # Jakes model but for scaterrers
            MAX_DOPPLER = max(shifts_valid)
            MIN_DOPPLER = min(shifts_valid)
            DELTA = (MAX_DOPPLER - MIN_DOPPLER) / n_bins

            bins = np.zeros(n_bins, dtype=np.complex128)
            bins_labels = np.linspace(MIN_DOPPLER, MAX_DOPPLER, num=n_bins + 1) + DELTA / 2
            bins_labels = bins_labels[:-1].copy()

            for _ in range(iterations):
                angle_perp, angle_tang = random_angles(random_type, ANGULAR_WIDTH_OF_PATCH)
                # Check if within circle
                if np.sqrt(angle_tang ** 2 + angle_perp ** 2) < ANGULAR_WIDTH_OF_PATCH / 2:
                    #  Limit receiver angle
                    if cos_of_angle_between(unit_vector_at_angle(RX_angle), [angle_tang, angle_perp]) > receiver_angle(
                            RX_width):
                        angle_perp += SATELLITE_ANGLE_PERPENDICULAR_FOR_CENTRE
                        angle_tang += SATELLITE_ANGLE_TANGENTIAL_FOR_CENTRE
                        doppler = doppler_shift(angle_perp, angle_tang)
                        bin_target = int((doppler - MIN_DOPPLER) // DELTA)
                        if bin_target >= n_bins:
                            bin_target = n_bins - 1
                        if bin_target <= -1:
                            bin_target = 0
                        if RAYLEIGH:
                            bins[bin_target] += random_unit_circle_complex_value()/total_signal_delay(angle_perp, angle_tang)/total_signal_delay(angle_perp, angle_tang)*4
                        else:
                            bins[bin_target] += 1/total_signal_delay(angle_perp, angle_tang)/total_signal_delay(angle_perp, angle_tang)*4
            if RX_angle in RX_angles_export and RX_width in RX_widths_export:
                file = f"model"
                np.savez(f"{main_folder}/{folder}/{subfolder}/{file}", bins=bins, bins_labels=bins_labels)

            bins = np.abs(bins)
            if not RAYLEIGH:
                bins = np.sqrt(bins)

            spreads.append(spread_from_bins(bins_labels, bins))

            # Plot the combined figure
            plt.figure(0)
            plt.plot(bins_labels, bins)

            # Plot the bar plots and save to folder
            plt.figure(random.randint(5, 10000000000))
            plt.bar(bins_labels, bins, width=DELTA)
            plt.xlabel("Hz")
            plt.ylabel("Amplitude")
            plt.title(f"Doppler spectrum for RXwidth={RX_width} and RXangle={int(RX_angle)}")
            plt.title(f"Doppler spectrum")

            file = f"RXangle_{int(RX_angle)}_Width_{int(RX_width)}.png"
            plt.savefig(f"{main_folder}/{folder}/{subfolder}/{file}")
            # plt.show(block=False)
            plt.close()
        plt.figure(3)
        plt.xlabel("RX width[degrees]")
        plt.ylabel("Doppler spread[Hz]")
        plt.title(f"Doppler spreads for RXangle = {RX_angle}")
        plt.title(f"Doppler spread vs beam width")
        plt.tight_layout()
        plt.plot(RX_widths, spreads)
        plt.savefig(f"{main_folder}/{folder}/{subfolder}/RXangle_{int(RX_angle)}_spreads.png", bbox_inches='tight')
        plt.clf()

        plt.figure(4)
        plt.xlabel("Number of antenna")
        plt.ylabel("Doppler spread[Hz]")
        plt.title(f"Doppler spread vs number of antenna")
        plt.plot(1/((np.sin(np.array(RX_widths)/4*2*np.pi/360))**2), spreads)
        plt.xlim(0, 270)
        plt.tight_layout()
        plt.savefig(f"{main_folder}/{folder}/{subfolder}/RXangle_{int(RX_angle)}_antenna.png", bbox_inches='tight')
        plt.clf()

        plt.figure(0)
        plt.savefig(f"{main_folder}/{folder}/{subfolder}/RXangle_{int(RX_angle)}.png")
        plt.clf()
        print(f"Successfully saved to: {main_folder}/" + folder + "/" + subfolder)
        #plt.show()

