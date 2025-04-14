import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import GroupContributionMethod as gcm

# -----------------------------------------------------------------------------
# Calculate mixture properties from the group contribution properties
# -----------------------------------------------------------------------------

# Fuel for GCM and data for validation (see fuelData/propertiesData for options)
# Options: 'decane', 'dodecane', 'heptane', 'posf10264', 'posf10289', 'posf10325'
fuel_name = 'posf10264'

# How much fuel? Amount does not matter.
drop_volume = 1e-9 # assume 1 microliter of fuel

# Plotting controls
fsize = 14 
ticksize = 12
line_thickness = 2
marker_size = 50

def drop_mass(fuel, Yi, T):
    # Determine the mass in kg of each component of the fuel
    return drop_volume / (fuel.molar_liquid_vol(T) @ Yi) * Yi * fuel.MW  # (kg)

# Get the fuel properties based on the GCM
fuel = gcm.groupContribution(fuel_name)

# Initial liquid mass fractions
Y_li = fuel.Y_0

# Mapping for fuel name to data and legend name
fuel_to_data = {
    'heptane': ('heptane-NIST.csv', 'NIST Data'),
    'decane':  ('decane-NIST.csv', 'NIST Data'),
    'dodecane': ('dodecane-NIST.csv', 'NIST Data'),
    'posf10264': ('posf10264.csv', 'Edwards et al. AFRL'),
    'posf10289': ('posf10289.csv', 'Edwards et al. AFRL'),
    'posf10325': ('posf10325.csv', 'Edwards et al. AFRL')
}
data_file, data_source = fuel_to_data.get(fuel_name)

# Load the experimental data
dataPath = os.path.join(fuel.fuelDataDir, "propertiesData")
data = pd.read_csv(os.path.join(dataPath, data_file), skiprows=[1])

# Separate properties and temperatures from data
T_nu_data = data.Temperature[data.Viscosity.notna()]
nu_data = data.Viscosity.dropna()
T_rho_data = data.Temperature[data.Density.notna()]
rho_data = data.Density.dropna()
T_pv_data = data.Temperature[data.VaporPressure.notna()]
pv_data = data.VaporPressure.dropna()
T_st_data = data.Temperature[data.SurfaceTension.notna()]
st_data = data.SurfaceTension.dropna() * 1e3
if 'posf' not in fuel_name:
    T_tc_data = data.Temperature[data.ThermalCond.notna()]
    tc_data = data.ThermalCond.dropna() 

# Generate temperature vectors and initialize GCM property vectors
T_rho = gcm.C2K(np.linspace(min(T_rho_data), max(T_rho_data), 100))
rho = np.zeros_like(T_rho)
T_nu = gcm.C2K(np.linspace(min(T_nu_data), max(T_nu_data), 100))
nu = np.zeros_like(T_nu)
T_pv = gcm.C2K(np.linspace(min(T_pv_data), max(T_pv_data), 100))
pv = np.zeros_like(T_pv)
T_st = gcm.C2K(np.linspace(min(T_st_data), max(T_st_data), 100))
st = np.zeros_like(T_st)
if 'posf' not in fuel_name:
    T_tc = gcm.C2K(np.linspace(min(T_tc_data), max(T_tc_data), 100))
    tc = np.zeros_like(T_tc)

# Calculate GCM properties for a range of temperatures
for i in range(len(T_rho)): 
    rho[i] = fuel.mixture_density(fuel.Y_0, T_rho[i]) # kg/m^3
    rho[i] *= 1e-3 # Convert to g/cm^3
for i in range(len(T_nu)): 
    mass = drop_mass(fuel, Y_li, T_nu[i])
    nu[i] = fuel.mixture_kinematic_viscosity(mass, T_nu[i]) # m^2/s
    nu[i] *= 1e6  # Convert to mm^2/s
for i in range(len(T_pv)): 
    mass = drop_mass(fuel, Y_li, T_pv[i])
    pv[i] = fuel.mixture_vapor_pressure(mass, T_pv[i]) # Pa
    pv[i] *= 1e-3  # Convert to kPa
for i in range(len(T_st)): 
    mass = drop_mass(fuel, Y_li, T_st[i])
    st[i] = fuel.mixture_surface_tension(mass,T_st[i],"Brock-Bird") # N/m
    st[i] *= 1e3 # Convert to dyne/cm
if 'posf' not in fuel_name:
    for i in range(len(T_tc)): 
        tc[i] = fuel.thermal_conductivity(T_tc[i])

# Create subplots
fig, axes = plt.subplots(1,4, figsize=(20, 5))
if "posf" in fuel_name:
    title_name = fuel_name.upper()
else:
    title_name = fuel_name.capitalize()

fig.suptitle(f"Fuel: {title_name}", fontsize=fsize + 2)

# Density
axes[0].plot(gcm.K2C(T_rho), rho, '-', label='Model Prediction', linewidth=line_thickness)
axes[0].scatter(T_rho_data, rho_data, label=data_source, facecolors='black', s=marker_size)
axes[0].set_xlim([min(T_rho_data-10), max(T_rho_data+10)])
axes[0].set_xlabel('Temperature (°C)', fontsize=fsize)
axes[0].set_ylabel('Density (g/cm³)', fontsize=fsize)
axes[0].tick_params(axis='both', labelsize=ticksize)

# Vapor Pressure
axes[1].plot(gcm.K2C(T_pv), pv, '-', label='Model Prediction', linewidth=line_thickness)
axes[1].scatter(T_pv_data, pv_data, label=data_source, facecolors='black', s=marker_size)
axes[1].set_xlim([min(T_pv_data-10), max(T_pv_data+10)])
axes[1].set_ylim([-1, max(pv_data)+ 1])
yticks = axes[1].get_yticks()
new_yticks = [tick for tick in yticks if tick >= 0]
axes[1].set_yticks(new_yticks)
axes[1].set_xlabel('Temperature (°C)', fontsize=fsize)
axes[1].set_ylabel('Vapor Pressure (kPa)', fontsize=fsize)
axes[1].tick_params(axis='both', labelsize=ticksize)


# Viscosity
axes[2].plot(gcm.K2C(T_nu), nu, '-', label='Model Prediction', linewidth=line_thickness)
axes[2].scatter(T_nu_data, nu_data, label=data_source, facecolors='black', s=marker_size)
axes[2].set_xlim([min(T_nu_data-10), max(T_nu_data+10)])
axes[2].set_xlabel('Temperature (°C)', fontsize=fsize)
axes[2].set_ylabel('Viscosity (mm²/s)', fontsize=fsize)
axes[2].tick_params(axis='both', labelsize=ticksize)

# Surface Tension
axes[3].plot(gcm.K2C(T_st), st, '-', label='Model Prediction', linewidth=line_thickness)
axes[3].scatter(T_st_data, st_data, label=data_source, facecolors='black', s=marker_size)
axes[3].set_xlim([min(T_st_data-10), max(T_st_data+10)])
axes[3].set_xlabel('Temperature (°C)', fontsize=fsize)
axes[3].set_ylabel('Surface Tension (dyne/cm)', fontsize=fsize)
axes[3].tick_params(axis='both', labelsize=ticksize)
axes[3].legend(fontsize=fsize-2)


plt.tight_layout(rect=[0.001, 0.03, 1, 0.99])
plt.show()

if 'posf' not in fuel_name:
    fig2 = plt.figure(figsize=(6, 5))
# Thermal Conductivity
    plt.plot(gcm.K2C(T_tc), tc, '-', label='Model Prediction', linewidth=line_thickness)
    plt.scatter(T_tc_data, tc_data, label=data_source, facecolors='black', s=marker_size)
    plt.xlim([min(T_tc_data-10), max(T_tc_data+10)])
    plt.xlabel('Temperature (°C)', fontsize=fsize)
    plt.ylabel('Thermal Conductivity (W/m/K)', fontsize=fsize)
    plt.legend(fontsize=10)
    plt.tick_params(axis='both', labelsize=ticksize)
    plt.tight_layout()
    plt.show()