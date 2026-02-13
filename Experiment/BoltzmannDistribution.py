import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Boltzmann distribution for TWO conformers
# ============================================================
# We consider two conformers:
#   - Conformer 1 (low energy): E1 = 0
#   - Conformer 2 (high energy): E2 = ΔE
#
# The Boltzmann formula gives:
#   P2/P1 = exp(-(E2-E1)/RT) = exp(-ΔE/RT)
#
# From this ratio we compute normalized populations:
#   P1 = 1/(1 + exp(-ΔE/RT))
#   P2 = 1 - P1
# ============================================================

# Temperature: 25 °C converted into Kelvin
T = 25 + 273.15

# Constant linking energy and temperature (molar statistical form)
R = 8.314462618  # J/(mol*K)

# Energy difference range: ΔE = E_high - E_low (kJ/mol)
dE_kJ = np.linspace(0, 10, 401)
dE_J = dE_kJ * 1000  # convert into Joules

# ------------------------------------------------------------
# 1) Population ratio P2/P1
# ------------------------------------------------------------
ratio = np.exp(-dE_J / (R * T))

# ------------------------------------------------------------
# 2) Normalized populations P1 and P2
# ------------------------------------------------------------
P1 = 1 / (1 + ratio)   # lower-energy conformer
P2 = 1 - P1            # higher-energy conformer

# ============================================================
# Plot A: Ratio P2/P1 (normal scale)
# ============================================================
plt.figure(figsize=(7.5,4.6))
plt.plot(dE_kJ, ratio)
plt.xlabel("Energy Difference ΔE (kJ/mol) = E_high - E_low")
plt.ylabel("Population ratio P2/P1")
plt.title("Boltzmann ratio for two conformers at 25 °C")
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# ============================================================
# Plot B: Ratio P2/P1 (log scale)
# Useful when ratios span several orders of magnitude
# ============================================================
plt.figure(figsize=(7.5,4.6))
plt.plot(dE_kJ, ratio)
plt.yscale("log")
plt.xlabel("Energy Difference ΔE (kJ/mol) = E_high - E_low")
plt.ylabel("Population ratio P2/P1 (log scale)")
plt.title("Boltzmann ratio (log scale) at 25 °C")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# ============================================================
# Plot C: Normalized populations P1 and P2
# ============================================================
plt.figure(figsize=(7.5,4.6))
plt.plot(dE_kJ, P1, label="P1 (low energy)")
plt.plot(dE_kJ, P2, label="P2 (high energy)")
plt.xlabel("Energy Difference ΔE (kJ/mol)")
plt.ylabel("Normalized population")
plt.title("Normalized conformer populations at 25 °C")
plt.grid(True, linestyle="--", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# Plot D: Textbook-style curve (% of conformers)
# With Y-axis labels: 50/50, 60/40, … 100/0
# ============================================================

# Convert population of the lower-energy conformer into %
P1_percent = 100 * P1

# Find ΔE where population is ~90% (90/10 split)
target = 90
idx = np.argmin(np.abs(P1_percent - target))
dE_90 = dE_kJ[idx]
P1_90 = P1_percent[idx]

plt.figure(figsize=(7.5,4.8))
plt.plot(dE_kJ, P1_percent)

# Custom Y ticks like the reference image
yticks = [50,60,70,80,90,100]
yticklabels = ["50/50","60/40","70/30","80/20","90/10","100/0"]
plt.yticks(yticks, yticklabels)

plt.xlabel("Energy Difference (kJ/mol)")
plt.ylabel("% of Conformers")
plt.title("Boltzmann population for 2 conformers at 25 °C")

# Dashed guide lines for 90/10
plt.axhline(target, linestyle="--")
plt.axvline(dE_90, linestyle="--")
plt.scatter([dE_90], [P1_90])

plt.ylim(50, 100)
plt.xlim(0, 10)
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

print("ΔE for ~90/10 population split =", round(dE_90,2), "kJ/mol")
