import numpy as np
import matplotlib.pyplot as plt

# ============================================================
#  Dimostrazione dell'effetto di V nella funzione torsionale
# ============================================================
# In un force field, l'energia torsionale descrive la variazione
# energetica associata alla rotazione attorno a un legame singolo.
#
# Una forma semplificata è:
#
#   E(θ) = V * [ 1 + cos(nθ) ]
#
# dove:
#   θ = angolo torsionale (0 → 360°)
#   V = ampiezza della barriera torsionale (kJ/mol)
#   n = periodicità (quante ripetizioni in 360°)
#
# ============================================================


# ------------------------------------------------------------
# 1) Creiamo un insieme di valori dell'angolo torsionale θ
# ------------------------------------------------------------
# θ varia da 0 a 2π radianti (equivalente a 0–360°)
theta = np.linspace(0, 2*np.pi, 400)

# ------------------------------------------------------------
# 2) Scegliamo la periodicità n
# ------------------------------------------------------------
# n indica quante volte il profilo energetico si ripete
# durante una rotazione completa.
#
# Per un legame C–C tipico spesso n = 3:
# → tre minimi e tre massimi in 360°
n = 3


# ------------------------------------------------------------
# 3) Scegliamo diversi valori di V
# ------------------------------------------------------------
# V determina quanto alta è la barriera torsionale:
#   V piccolo → rotazione facile
#   V grande  → rotazione sfavorita
#
# Qui confrontiamo tre casi:
V_values = [1,0.3]  # kJ/mol
T_values = [0,40]  #

# ------------------------------------------------------------
# 4) Iniziamo la figura del grafico
# ------------------------------------------------------------
plt.figure(figsize=(7.5, 4.8))

# Convertiamo θ in gradi per renderlo più leggibile
theta_deg = theta * 180 / np.pi
# ------------------------------------------------------------
# 5) Calcoliamo e plottiamo E(θ) per ogni valore di V
# ------------------------------------------------------------
e=[]
for V, T in zip(V_values, T_values):

    # Energia torsionale secondo la formula:
    # E(θ) = V * (1 + cos(nθ))
    #
    # Il termine (1 + cos(...)) oscilla tra:
    #   0  (minimo energetico)
    #   2  (massimo energetico)
    #
    # Quindi l'energia varia tra:
    #   0  e 2V
    E = V * (1 + np.cos(n * theta + T))
    e.append(E)
    # Disegniamo la curva
    plt.plot(theta_deg, E, label=f"V = {V} kJ/mol")

# ------------------------------------------------------------
# 6) Aggiungiamo etichette e stile
# ------------------------------------------------------------
plt.xlabel("Torsional angle θ (degrees)")
plt.ylabel("Torsional energy E(θ)")
plt.title("Effect of V on torsional barrier height (n = 3)")
plt.grid(True, linestyle="--", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
