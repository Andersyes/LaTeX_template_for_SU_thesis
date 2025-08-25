import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("homogeneous_freezing_CNT.csv")

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(111)

# Plot J and rate per droplet vs temperature
ax1.set_yscale('log')
ax1.plot(df['Temperature_C'], df['J_m^-3_s^-1'], label='J (m$^{-3}$ s$^{-1}$)', color = 'firebrick')
ax1.set_xlabel('Temperature (°C)')
ax1.set_ylabel('Homogeneous nucleation rate J (m$^{-3}$ s$^{-1}$)')
ax1.grid(True, which='both', linestyle=':', linewidth=0.5)

# Vertical line at -38 °C
ax1.axvline(-38.0, linestyle='--', color='gray', linewidth=1.0)
ax1.text(-37.5, 1e38, '-38 °C', verticalalignment='top')

## Twin y-axis for freezing rate per droplet
#ax2 = ax1.twinx()
#ax2.set_yscale('log')
#ax2.plot(df['Temperature_C'], df['Rate_per_drop_s^-1'], label='Rate per drop (s$^{-1}$)', color='orange')
#ax2.set_ylabel('Freezing rate per droplet (s$^{-1}$)')

# Twin y-axis for mean wait
ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.plot(df['Temperature_C'], df['Mean_wait_s'], label='Mean wait before freezing', color='forestgreen')
ax2.set_ylabel('Time (s)')

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='lower left')

plt.title('Classical homogeneous freezing rates')
plt.tight_layout()
plt.savefig('./images/homogeneous_nucleation_rate.png', dpi=300)
plt.show()

