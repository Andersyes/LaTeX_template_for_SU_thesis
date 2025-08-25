# Plot full Delta G(r) vs r for liquid->ice 
import numpy as np
import matplotlib.pyplot as plt

# Physical constants and parameters
sigma = 0.030       # N/m (J/m^2) ice-water interfacial tension
rho_i = 917.0       # kg/m^3
Lf = 3.34e5         # J/kg
Tm = 273.15         # K
kB = 1.380649e-23   # J/K

# Temperature
T_C = -20.0
T = T_C + 273.15
DeltaT = Tm - T

# volumetric Gibbs free energy difference (linear approx)
Delta_g_v = - rho_i * Lf * (DeltaT / Tm)  # negative (J/m^3)
abs_Delta_g_v = abs(Delta_g_v)

# critical radius and barrier (homogeneous)
r_c = 2.0 * sigma / abs_Delta_g_v
DeltaG_c = (16.0 * np.pi / 3.0) * (sigma**3) / (abs_Delta_g_v**2)

# radius array: cover from near 0 to 5*r_c
r = np.linspace(1e-11, 5*r_c, 1000)

# Homogeneous Delta G(r)
DeltaG_r_hom = 4.0 * np.pi * r**2 * sigma + (4.0/3.0) * np.pi * r**3 * Delta_g_v

# Heterogeneous geometric factor (spherical cap approximation)
def f_theta(theta_rad):
    return ((1 - np.cos(theta_rad))**2 * (2 + np.cos(theta_rad))) / 4.0

theta_degs = [30, 60, 90, 120]   # chosen contact angles
theta_rads = np.radians(theta_degs)
f_vals = [f_theta(th) for th in theta_rads]

# Approximate heterogeneous DeltaG(r) by scaling the homogeneous DeltaG(r) by f(theta)
# Note: this is an approximation that correctly scales the barrier height; exact shape differs.
DeltaG_r_het = {deg: f * DeltaG_r_hom for deg, f in zip(theta_degs, f_vals)}

## Plot DeltaG(r) in Joules (linear y-scale to show hump; we can also include kT scale)
#plt.figure(figsize=(8,5))
#plt.plot(r*1e9, DeltaG_r_hom, label=f'Homogeneous (r_c={r_c*1e9:.2f} nm)')
#for deg in theta_degs:
#    plt.plot(r*1e9, DeltaG_r_het[deg], label=f'Hetero θ={deg}° (f={f_theta(np.radians(deg)):.3f})')
#
## Mark critical radius
#plt.axvline(r_c*1e9, color='k', linestyle='--', linewidth=1)
#plt.text(r_c*1e9*1.05, max(DeltaG_r_hom)*0.6, f'r_c = {r_c*1e9:.2f} nm', rotation=90, verticalalignment='center')
#
#plt.xlabel('Embryo radius r (nm)')
#plt.ylabel('ΔG(r) (J)')
#plt.title(f'Gibbs free energy ΔG(r) vs embryo radius at {T_C:.0f} °C\n(liquid → ice, CNT; hetero approximated by f(θ) scaling)')
#plt.legend()
#plt.grid(True, linestyle=':', linewidth=0.5)
#plt.ylim(bottom=0)
#plt.xlim(left=0, right=5*r_c*1e9)
#plt.tight_layout()
#plt.show()

colors = ['olivedrab','mediumaquamarine','cadetblue','teal']
lstyles= ['-','--','-.',':']
# Also plot in units of kB T for clarity (log scale for y to span range)
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.plot(r*1e9, DeltaG_r_hom/(kB*T), linewidth=2, color = 'black', label='Homogeneous')
i = -1
for deg in theta_degs:
    i+=1
    ax.plot(r*1e9, DeltaG_r_het[deg]/(kB*T), linestyle = lstyles[i], color=colors[i], label=f'Hetero θ={deg}° (f={f_theta(np.radians(deg)):.3f})')

# plot r_c
ax.vlines(r_c*1e9,0,DeltaG_c/(kB*T), color='k', linestyle=':', linewidth=1)
ax.text(r_c*1e9*1.03, 270, r'$r_c = {}$ nm'.format(round(r_c*1e9,2)), rotation=0, fontsize = 12, verticalalignment='center')
#ax.hlines(DeltaG_c/(kB*T), 0, r_c*1e9, color='k', linestyle=':', linewidth=1)



ax.set_xlabel('Embryo radius r (nm)', fontsize = 15)
ax.set_ylabel(r'$\Delta G(r) / (k_B T)$', fontsize = 15)
#ax.set_title(r'$\Delta G(r) / (k_B T)$ vs embryo radius at {}°C'.format(T_C), fontsize = 13)
ax.set_title(r'Energy barrier for ice nucleation from liquid water', fontsize = 15)
ax.legend(fontsize = 12)
ax.grid(True, which='both', linestyle=':', linewidth=0.5)
ax.tick_params(which='both', labelsize = 12)
ax.set_xlim(left=0, right=4.2)
ax.set_ylim(1e-2, 300)
plt.tight_layout()
plt.savefig('./images/nucleation_barrier.png', dpi=300)

# Print key numbers
print(f"T = {T_C} °C (T = {T:.2f} K)")
print(f"DeltaT = {DeltaT:.2f} K")
print(f"Critical radius r_c = {r_c:.3e} m = {r_c*1e9:.2f} nm")
print(f"Homogeneous barrier DeltaG_c = {DeltaG_c:.3e} J = {DeltaG_c/(kB*T):.3e} k_B T")
