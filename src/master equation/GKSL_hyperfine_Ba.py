# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 13:56:15 2026

@author: USER-01
"""

import numpy as np
import library
from qutip import *
from scipy.constants import m_p, c, h, hbar, epsilon_0
import matplotlib.pyplot as plt

# ── Fine-structure states ────────────────────────────────────────────────────
state_1 = library.State(configuration="6s2",  S=0, L=0, J=0, energy=0)
state_2 = library.State(configuration="6s6p", S=0, L=1, J=1, energy=18060.261)
state_3 = library.State(configuration="6s6d", S=0, L=2, J=2, energy=30236.826)
J_states = [state_1, state_2, state_3]

# ── Hyperfine states (I = 3/2) ───────────────────────────────────────────────
# state_1 (J=0): F = 3/2          → 1 state   [index 0]
# state_2 (J=1): F = 1/2, 3/2, 5/2 → 3 states [indices 1-3]
# state_3 (J=2): F = 1/2, 3/2, 5/2, 7/2 → 4 states [indices 4-7]
I = 3/2
A_hf = [-109.5, -109.5, -109.5]   # placeholder A coefficients (MHz) for each J level
B_hf = [0.0,    50.09,  50.09]    # placeholder B coefficients

hyperfine_states = []
for idx, state in enumerate(J_states):
    hf = state.hyperfine_states(I=I, A=A_hf[idx], B=B_hf[idx])
    hyperfine_states.append(hf)

# Flatten into a single ordered list: [g0 | e1,e2,e3 | d1,d2,d3,d4]
all_states = hyperfine_states[0] + hyperfine_states[1] + hyperfine_states[2]
N = len(all_states)
print(N)

# Map state object → index
state_index = {id(s): i for i, s in enumerate(all_states)}

# ── Fine-structure transitions (for laser frequencies & A coefficients) ──────
transition_12 = library.Transition(state_1, state_2, linewidth=20e6, A=1.19e8)
transition_23 = library.Transition(state_2, state_3, linewidth=20e6, A=5.00e6)
transition_13 = library.Transition(state_1, state_3, linewidth=20e6, A=1e5)

laser_1 = library.Laser(frequency=transition_12.frequency(), intensity=1e5,  linewidth=4e6)
laser_2 = library.Laser(frequency=transition_23.frequency(), intensity=1e5,  linewidth=4e6)
laser_3 = library.Laser(frequency=c / 405e-9,               intensity=1e10, linewidth=4e6)

# Laser frequencies (angular, rad/s)
nu1 = laser_1.frequency * 2 * np.pi
nu2 = laser_2.frequency * 2 * np.pi
nu3 = laser_3.frequency * 2 * np.pi   # ionisation laser – not included in H

# ── Helper: assign laser driving a pair of J-manifolds ───────────────────────
def laser_for_pair(i_group, j_group):
    """Return (laser, nu_laser) for a given pair of J-manifold indices."""
    if (i_group, j_group) in [(0, 1), (1, 0)]:
        return laser_1, nu1
    elif (i_group, j_group) in [(1, 2), (2, 1)]:
        return laser_2, nu2
    else:
        return None, None

# Group index for each state in all_states
def group_of(state_idx):
    if state_idx < len(hyperfine_states[0]):
        return 0
    elif state_idx < len(hyperfine_states[0]) + len(hyperfine_states[1]):
        return 1
    else:
        return 2

# ── Allowed hyperfine transitions (ΔF = 0, ±1) ──────────────────────────────
# Collect pairs (lower_idx, upper_idx, Omega_ij)
allowed_transitions = []

for i, si in enumerate(all_states):
    gi = group_of(i)
    
    for j, sj in enumerate(all_states):
        
        gj = group_of(j)
        
        if gj <= gi:
            continue   # only upper triangle (j > i energetically)
        if abs(gi - gj) != 1:
            continue   # only adjacent J-manifolds for dipole transitions
        
        # ΔF selection rule
        Fi = getattr(si, 'F', si.J)   # fall back to J if F not set
        Fj = getattr(sj, 'F', sj.J)
        if abs(Fj - Fi) > 1:
            continue

        laser, nu_laser = laser_for_pair(gi, gj)
        if laser is None:
            continue

        # Determine J and A coefficient of the fine-structure transition
        if (gi, gj) == (0, 1):
            A_coeff = transition_12.A
            J_lower, J_upper = state_1.J, state_2.J
        else:
            A_coeff = transition_23.A
            J_lower, J_upper = state_2.J, state_3.J

        # Transition dipole moment from A coefficient (fine-structure line)
        # d = sqrt(3 π ε0 ħ c³ A / ω³)
        omega_0 = (sj.energy - si.energy) * c * 100 * 2 * np.pi   # rad/s
        if omega_0 <= 0:
            continue
        d_ij = np.sqrt(3 * np.pi * epsilon_0 * hbar * c**3 * A_coeff / omega_0**3)

        # Clebsch-Gordan scaling
        # Use qutip clebsch(j1, j2, j3, m1, m2, m3)
        # For a π-transition (Δm_F = 0) as a placeholder
        # CG = <Fi, 0; 1, 0 | Fj, 0>  (m_F = 0 projection)
        cg = clebsch(Fi, 1, Fj, 0, 0, 0)

        E0 = laser.electric_field_amplitude()
        Omega_ij = d_ij * E0 * abs(cg) / hbar   # rad/s
        
        # print(Omega_ij)

        allowed_transitions.append((i, j, Omega_ij, nu_laser))

print(f"Number of allowed hyperfine transitions: {len(allowed_transitions)}")

# ── Build Hamiltonian in rotating frame ──────────────────────────────────────
# H = Σ_i ħ(ω_i − ν_i)|i><i|  +  (ħ/2) Σ_{(i,j)} (Ω_ij|i><j| + h.c.)
# ω_i: state frequency (rad/s);  ν_i: laser frequency driving state i

# State frequencies (rad/s) relative to ground
omega = np.array([s.energy * c * 100 * 2 * np.pi for s in all_states])

# Assign a laser frequency ν_i to each state for the diagonal
def laser_nu_for_state(idx):
    g = group_of(idx)
    if g == 0:
        return 0.0   # ground – reference
    elif g == 1:
        return nu1
    else:
        return nu1 + nu2   # two-photon accumulated detuning

nu_state = np.array([laser_nu_for_state(i) for i in range(N)])

H = Qobj(np.zeros((N, N), dtype=complex))

# Diagonal terms
for i in range(N):
    H += hbar * (omega[i] - nu_state[i]) * basis(N, i) * basis(N, i).dag()

# Off-diagonal laser interaction terms
for (i, j, Omega_ij, _) in allowed_transitions:
    # |i><j| + |j><i|  (i < j, i is lower)
    H += (hbar / 2) * Omega_ij * (basis(N, i) * basis(N, j).dag()
                                 + basis(N, j) * basis(N, i).dag())

print(H)


# ── Jump operators (Lindblad) ─────────────────────────────────────────────────
c_ops = []

for (i, j, _, _) in allowed_transitions:
    # Determine A coefficient for this pair
    gi, gj = group_of(i), group_of(j)
    
    if (gi, gj) == (0, 1):
        A_decay = transition_12.A
        
    else:
        A_decay = transition_23.A

    # Spontaneous emission: |i><j| (decay from j → i)
    decay_op = np.sqrt(A_decay) * basis(N, i) * basis(N, j).dag()
    c_ops.append(decay_op)

'''
# Ionisation from state_3 manifold (laser_3)
sigma_ion = 6e-21   # m²  (placeholder)
W_ion = sigma_ion * laser_3.intensity / (h * laser_3.frequency)

for idx in range(len(hyperfine_states[0]) + len(hyperfine_states[1]), N):
    ion_op = np.sqrt(W_ion) * basis(N, idx) * basis(N, idx).dag()
    c_ops.append(ion_op)
    
'''

# ── Initial state: ground hyperfine level ─────────────────────────────────────
psi0 = basis(N, 0)
rho0 = ket2dm(psi0)

# ── Time evolution ────────────────────────────────────────────────────────────
t_end = 1e-6        # 1 µs
n_steps = 500
tlist = np.linspace(0, t_end, num=n_steps)

result = mesolve(H, 
                 rho0, 
                 tlist, 
                 c_ops, 
                 e_ops=[basis(N, i) * basis(N, i).dag() for i in range(N)])

# ── Plot populations ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

group_labels = [
    (range(len(hyperfine_states[0])),                                                          "Ground (J=0) hyperfine"),
    (range(len(hyperfine_states[0]), len(hyperfine_states[0])+len(hyperfine_states[1])),       "Excited (J=1) hyperfine"),
    (range(len(hyperfine_states[0])+len(hyperfine_states[1]), N),                              "Upper (J=2) hyperfine"),
]

colors = plt.cm.tab10.colors

for ax, (indices, title) in zip(axes, group_labels):
    for idx in indices:
        F_val = getattr(all_states[idx], 'F', '?')
        ax.plot(tlist * 1e6, result.expect[idx], label=f"F={F_val}", color=colors[idx % 10])
        # print(result.expect[idx])
    ax.set_ylabel("Population")
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.set_ylim(-0.05, 1.05)

axes[-1].set_xlabel("Time (µs)")
plt.tight_layout()
plt.show()

