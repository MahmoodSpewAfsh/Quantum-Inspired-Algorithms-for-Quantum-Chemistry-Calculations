# Define a colormap for better color differentiation
colors = ['purple', 'blue', 'green', 'cyan', 'red']
markers = ['o', 's', 'D', 'x', '^']
linestyles = ['dashed', 'solid', 'dashdot', 'dotted', (0, (3, 5, 1, 5))]  # Added more variations

# Data series
energy_values = [energy2, energy3, energy4, energy5, energy1]
labels = ["r=2", "r=3", "r=4", "r=5", "r=6"]

fig, ax = plt.subplots(figsize=(12, 7))

# Plot each data series
for i, (energy, label) in enumerate(zip(energy_values, labels)):
    ax.plot(bond_lengths, energy, color=colors[i], linewidth=1.5,
            label=label, marker=markers[i], linestyle=linestyles[i])

# Plot exact energies
ax.plot(
    bond_lengths,
    exact_energies,
    color='black',
    linewidth=3,
    linestyle='solid',
    marker=None,
    label='Exact (CASCI(2,2))'
)
# Plot HF energies
ax.plot(
    bond_lengths,
    hf_energies,
    color='orange',
    linewidth=2.5,
    linestyle='solid',
    label='HF'
)
# Exact dissociation limit
ax.axhline(
    y=E_diss,
    color='black',
    linestyle='--',
    linewidth=3,
    label='Exact dissociation limit'
)
# Labels and styling
ax.set_xlabel('Bond Length', fontsize=12)
ax.set_ylabel('Energy of Hydrogen', fontsize=12)
ax.set_title('Hydrogen Energy vs. Bond Length with SFC', fontsize=14)
ax.legend()
ax.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
#Saving the Result
np.save("SFCH2r=2", energy2)
np.save("SFCH2r=3", energy3)
np.save("SFCH2r=4", energy4)
np.save("SFCH2r=5", energy5)
np.save("SFCH2r=6", energy1)
