# Computing homogeneous nucleation rate J(T) using classical nucleation theory
# and plotting J (m^-3 s^-1) and the droplet freezing rate J*Vdrop (s^-1) vs temperature.
import numpy as np
import matplotlib.pyplot as plt

# Physical constants and parameters (typical values)
sigma_iw = 0.030        # N/m, ice-water interfacial tension
rho_i = 917.0           # kg/m^3, ice density
Lf = 3.34e5             # J/kg, latent heat of fusion
Tm = 273.15             # K, melting temperature
kB = 1.380649e-23       # J/K, Boltzmann constant
J0 = 1e35               # m^-3 s^-1, kinetic prefactor (uncertain)
r_drop = 10e-6          # m, droplet radius (10 micron)
V_drop = 4/3 * np.pi * r_drop**3  # m^3

# Temperature range (Kelvin) - from -55°C to -30°C
T_C = np.linspace(-55.0, -30.0, 501)
T = T_C + 273.15

# Helper: volumetric Gibbs free energy difference approx
# Delta g_v ≈ rho_i * Lf * (Tm - T) / Tm
DeltaT = Tm - T
Delta_g_v = rho_i * Lf * (DeltaT / Tm)  # J/m^3

# Avoid division by zero or tiny numbers for T >= Tm
Delta_g_v[Delta_g_v <= 0] = 1e-30

# Critical free energy barrier (classical formula)
DeltaGc = (16.0 * np.pi / 3.0) * (sigma_iw**3) / (Delta_g_v**2)  # J

# Homogeneous nucleation rate
exponent = -DeltaGc / (kB * T)
# To avoid underflow for extremely negative exponents, clip exponent at a reasonable min
exponent_clipped = np.clip(exponent, -700, 700)
J = J0 * np.exp(exponent_clipped)  # m^-3 s^-1

# Rate per droplet and mean waiting time (seconds)
rate_per_drop = J * V_drop  # s^-1
with np.errstate(divide='ignore', invalid='ignore'):
    mean_wait_s = np.where(rate_per_drop > 0, 1.0 / rate_per_drop, np.inf)

fig, ax1 = plt.subplots(figsize=(8,5))

# Plot: J and J*V_drop vs Temperature (°C)
#ax1.set_yscale('log')
#ax1.plot(T_C, J, label='J (m$^{-3}$ s$^{-1}$)', color='firebrick')
#ax1.set_xlabel('Temperature (°C)')
#ax1.set_ylabel('Homogeneous nucleation rate J (m$^{-3}$ s$^{-1}$)')
#ax1.set_ylim(1e-2, 1e40)
## Twin axis for rate per droplet
#ax2 = ax1.twinx()
#ax2.set_yscale('log')
#ax2.plot(T_C, rate_per_drop, label='J * V_drop (s$^{-1}$)', color = 'orange')
#ax2.set_ylabel('Freezing rate per droplet (s$^{-1}$)')
#ax2.set_ylim(1e-20, 1e5)
ax1.grid(True, which='both', linestyle=':', linewidth=0.5)

# mark -38°C
ax1.axvline(-38.0, linestyle='--', linewidth=1.0, color = 'black')
ax1.text(-37.5, 5e1, '-38 °C', verticalalignment='top', color = 'black', fontsize = 12)

ax1.set_yscale('log')
ax1.plot(T_C, rate_per_drop, label=r'J * $V_{drop}$ (s$^{-1}$)', color = 'darkgreen')
ax1.set_ylabel('Freezing rate per droplet (s$^{-1}$)', fontsize=13)
ax1.set_xlabel('Temperature (°C)', fontsize=13)
ax1.tick_params(which='both', labelsize = 11)
ax1.set_ylim(1e-23, 1e3)

# Twin y-axis for mean wait
ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.plot(T_C, mean_wait_s, label='Mean wait before freezing', color='darkcyan')
ax2.set_ylabel('Time (s)', fontsize=13)
ax2.set_ylim(5e-2, 1e29)

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='center left')

plt.title('Classical homogeneous freezing: droplet freezing rate and mean time before freezing\n(σ = 0.030 N/m, J0 = 1e35, r = 10 µm)', fontsize=13)
plt.tight_layout()
plt.savefig('./images/homogeneous_nucleation_rate.png', dpi=300)

# Also provide a small table of selected temperatures
import pandas as pd
sel_temps = np.array([-50, -45, -40, -39, -38, -37, -36])
rows = []
for tc in sel_temps:
    idx = np.argmin(np.abs(T_C - tc))
    rows.append({
        'T (°C)': tc,
        'J (m^-3 s^-1)': f"{J[idx]:.3e}",
        'rate per drop (s^-1)': f"{rate_per_drop[idx]:.3e}",
        'mean wait (s)': f"{mean_wait_s[idx]:.3e}" if np.isfinite(mean_wait_s[idx]) else "inf"
    })
df = pd.DataFrame(rows)
import caas_jupyter_tools as cjt; cjt.display_dataframe_to_user("Selected CNT values", df)

