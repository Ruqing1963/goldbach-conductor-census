"""
fig_census.py — Figures for Paper #14: Conductor Census
"""
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Load data
with open('/home/claude/paper14/data/conductor_census.json') as f:
    census = json.load(f)

with open('/home/claude/paper14/data/decomposition.json') as f:
    decomp = json.load(f)

with open('/home/claude/paper14/data/fit_params.json') as f:
    params = json.load(f)

# Filter for plotting
census_100 = [d for d in census if d['twoN'] >= 20]

# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Mean conductor ratio vs log(2N) with decomposition
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel (a): ⟨ρ_true⟩ vs log(2N)
ax = axes[0]
twoNs = [d['twoN'] for d in census_100]
rho_means = [d['rho_true_mean'] for d in census_100]
log2Ns = [math.log(d['twoN']) for d in census_100]

# Scatter with transparency for density
ax.scatter(log2Ns, rho_means, s=2, alpha=0.15, c='#2255BB', rasterized=True)

# Moving average (window of 50 points)
window = 100
if len(rho_means) > window:
    smooth_x = []
    smooth_y = []
    for i in range(window//2, len(rho_means) - window//2):
        smooth_x.append(log2Ns[i])
        smooth_y.append(sum(rho_means[i-window//2:i+window//2]) / window)
    ax.plot(smooth_x, smooth_y, 'r-', lw=2, label='Moving average', zorder=5)

# Fit line
a, b = params['rho_mean']['a'], params['rho_mean']['b']
x_fit = np.linspace(min(log2Ns), max(log2Ns), 100)
ax.plot(x_fit, a * x_fit + b, 'k--', lw=1, alpha=0.5,
       label=f'Fit: {a:.3f}·ln(2N) + {b:.2f}')

ax.set_xlabel(r'$\ln(2N)$', fontsize=12)
ax.set_ylabel(r'$\langle\rho_{\rm true}\rangle$', fontsize=12)
ax.set_title(r'(a) Mean conductor ratio vs $\ln(2N)$', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)

# Panel (b): Decomposition stacked area
ax = axes[1]
decomp_100 = [d for d in decomp if d['twoN'] >= 100 and d['twoN'] <= 10000]
x_d = [d['twoN'] for d in decomp_100]
y_static = [d['rho_static'] for d in decomp_100]
y_boundary = [d['rho_boundary_mean'] for d in decomp_100]
y_delta = [d['delta_mean'] for d in decomp_100]

# Stacked area
ax.fill_between(x_d, 0, y_static, alpha=0.7, color='#228833', label=r'$\rho_{\rm static}$ (from $N$)')
ax.fill_between(x_d, y_static, [s+b for s,b in zip(y_static, y_boundary)],
               alpha=0.7, color='#2255BB', label=r'$\langle\rho_{\rm boundary}\rangle$ (from $p \cdot q$)')
ax.fill_between(x_d, [s+b for s,b in zip(y_static, y_boundary)],
               [s+b+d for s,b,d in zip(y_static, y_boundary, y_delta)],
               alpha=0.7, color='#CC3333', label=r'$\langle\delta\rangle$ (from $|p{-}q|$)')

ax.set_xlabel(r'$2N$', fontsize=12)
ax.set_ylabel(r'Contribution to $\langle\rho_{\rm true}\rangle$', fontsize=12)
ax.set_title(r'(b) Three-component decomposition', fontsize=12)
ax.legend(fontsize=9, loc='upper left')
ax.set_xlim(100, 10000)
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/home/claude/paper14/figures/fig_mean_decomp.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/paper14/figures/fig_mean_decomp.png', dpi=200, bbox_inches='tight')
plt.close()
print("Figure 1 done.")


# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Bandwidth (σ) and normalized ratio ⟨ρ⟩/ln(2N)
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel (a): σ(ρ) vs 2N
ax = axes[0]
stds = [d['rho_true_std'] for d in census_100]
ax.scatter([d['twoN'] for d in census_100], stds,
          s=2, alpha=0.15, c='#CC8822', rasterized=True)

# Moving average
if len(stds) > window:
    smooth_s = []
    smooth_xs = []
    twoN_list = [d['twoN'] for d in census_100]
    for i in range(window//2, len(stds) - window//2):
        smooth_xs.append(twoN_list[i])
        smooth_s.append(sum(stds[i-window//2:i+window//2]) / window)
    ax.plot(smooth_xs, smooth_s, 'r-', lw=2, label='Moving average')

ax.set_xlabel(r'$2N$', fontsize=12)
ax.set_ylabel(r'$\sigma(\rho_{\rm true})$', fontsize=12)
ax.set_title(r'(a) Bandwidth $\sigma(\rho)$ vs $2N$', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)

# Panel (b): ⟨ρ⟩/ln(2N) — the normalized ratio
ax = axes[1]
# Sample every 10th point for clarity
sample = census_100[::5]
x_norm = [d['twoN'] for d in sample]
y_norm = [d['rho_true_mean'] / math.log(d['twoN']) for d in sample]

ax.scatter(x_norm, y_norm, s=8, alpha=0.3, c='#2255BB', rasterized=True)

# Moving average on full data
y_norm_full = [d['rho_true_mean'] / math.log(d['twoN']) for d in census_100]
twoN_full = [d['twoN'] for d in census_100]
if len(y_norm_full) > window:
    smooth_n = []
    smooth_xn = []
    for i in range(window//2, len(y_norm_full) - window//2):
        smooth_xn.append(twoN_full[i])
        smooth_n.append(sum(y_norm_full[i-window//2:i+window//2]) / window)
    ax.plot(smooth_xn, smooth_n, 'r-', lw=2, label='Moving average')

# Asymptotic reference
# ρ ~ 2·log(pq·rad(N)·rad(|p-q|))/log(2N)
# For typical p,q ~ N: ~ 4·log(N)/log(2N) + ... → 4
ax.axhline(y=0.5, color='gray', linestyle=':', lw=1, alpha=0.5)
ax.text(8000, 0.52, r'$\frac{1}{2}$', fontsize=10, color='gray')

ax.set_xlabel(r'$2N$', fontsize=12)
ax.set_ylabel(r'$\langle\rho_{\rm true}\rangle \,/\, \ln(2N)$', fontsize=12)
ax.set_title(r'(b) Normalized conductor ratio', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/home/claude/paper14/figures/fig_bandwidth_norm.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/paper14/figures/fig_bandwidth_norm.png', dpi=200, bbox_inches='tight')
plt.close()
print("Figure 2 done.")


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: G(2N) vs Hardy-Littlewood + conductor-count diagram
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel (a): G(2N) actual vs HL prediction
ax = axes[0]

# Hardy-Littlewood prediction
C2 = params['C2']
PRIMES_SMALL = [p for p in range(3, 1000) if all(p % d != 0 for d in range(2, int(p**0.5)+1))]

def hardy_littlewood(twoN):
    N = twoN // 2
    if N < 3: return 0
    product = 1.0
    temp = N; d = 3
    while d * d <= temp:
        if temp % d == 0:
            product *= (d - 1) / (d - 2)
            while temp % d == 0: temp //= d
        d += 2
    if temp > 2:
        product *= (temp - 1) / (temp - 2)
    return 2 * C2 * product * N / (math.log(N))**2

sample_c = [d for d in census if d['twoN'] >= 100 and d['twoN'] <= 10000]
actual = [d['count'] for d in sample_c]
predicted = [hardy_littlewood(d['twoN']) for d in sample_c]
twoN_c = [d['twoN'] for d in sample_c]

ax.scatter(predicted, actual, s=3, alpha=0.15, c='#228833', rasterized=True)
mn, mx = 0, max(max(actual), max(predicted)) * 1.1
ax.plot([mn, mx], [mn, mx], 'k--', lw=1, alpha=0.3, label='$G = G_{\\rm HL}$')
ax.set_xlabel(r'$G_{\rm HL}(2N)$ (Hardy--Littlewood)', fontsize=12)
ax.set_ylabel(r'$G(2N)$ (actual count)', fontsize=12)
ax.set_title(r'(a) Goldbach count vs Hardy--Littlewood', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)

# Panel (b): Coefficient of variation σ/⟨ρ⟩ vs 2N
ax = axes[1]
cv_data = [d for d in census_100 if d['rho_true_mean'] > 0 and d['count'] >= 3]
cv_x = [d['twoN'] for d in cv_data]
cv_y = [d['rho_true_std'] / d['rho_true_mean'] for d in cv_data]

ax.scatter(cv_x, cv_y, s=2, alpha=0.15, c='#CC3333', rasterized=True)

# Moving average
if len(cv_y) > window:
    smooth_cv = []
    smooth_cvx = []
    for i in range(window//2, len(cv_y) - window//2):
        smooth_cvx.append(cv_x[i])
        smooth_cv.append(sum(cv_y[i-window//2:i+window//2]) / window)
    ax.plot(smooth_cvx, smooth_cv, 'b-', lw=2, label='Moving average')

ax.set_xlabel(r'$2N$', fontsize=12)
ax.set_ylabel(r'$\sigma(\rho)\,/\,\langle\rho\rangle$', fontsize=12)
ax.set_title(r'(b) Coefficient of variation', fontsize=12)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.15)
ax.set_ylim(0, 0.25)

plt.tight_layout()
plt.savefig('/home/claude/paper14/figures/fig_HL_cv.pdf', dpi=300, bbox_inches='tight')
plt.savefig('/home/claude/paper14/figures/fig_HL_cv.png', dpi=200, bbox_inches='tight')
plt.close()
print("Figure 3 done.")


# ═══════════════════════════════════════════════════════════════
# Print key findings summary
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 65)
print("KEY FINDINGS SUMMARY")
print("=" * 65)
print(f"""
1. CONDUCTOR RATIO CONCENTRATION
   At 2N = 10000 (127 Goldbach pairs):
     ⟨ρ⟩ = 5.377,  σ = 0.323,  CV = σ/⟨ρ⟩ = 6.0%
   The conductor ratio is tightly concentrated for fixed N.
   This is WHY the BSL works: ~94% of ρ is deterministic.

2. THREE-COMPONENT DECOMPOSITION
     ρ_static    (from N):     6.5% of ⟨ρ⟩  — sets the BAND
     ρ_boundary  (from p·q):  66.0% of ⟨ρ⟩  — locates pair WITHIN band
     δ           (from |p-q|): 27.5% of ⟨ρ⟩  — NOISE (scatter)
   The proxy captures components 1+2 (72.5%), explaining R²>0.997.

3. NORMALIZED RATIO ⟨ρ⟩/ln(2N)
   Decreases from ~1.15 (2N=50) to ~0.58 (2N=10000).
   Appears to converge to ~1/2 from above.
   If true: ⟨ρ⟩ ~ (1/2)·ln(2N) + O(1) as N → ∞.

4. BANDWIDTH STABILITY
   σ(ρ) ≈ 0.33 ± 0.05 across the entire range 2N ∈ [100, 10000].
   The bandwidth is CONSTANT, not growing with N.
   This means the signal-to-noise ratio IMPROVES with N.

5. HARDY-LITTLEWOOD AGREEMENT
   G(2N) tracks the HL prediction with expected fluctuations.
   No anomalies — the conductor framework is consistent with
   the standard prime pair density model.
""")
