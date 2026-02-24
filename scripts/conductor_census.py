"""
conductor_census.py — Large-scale conductor ratio census for Goldbach-Frey curves
Direction A+D: Computation + Analytic Number Theory Bridge

For each even 2N, enumerate ALL Goldbach pairs (p,q), compute:
  - ρ_true = 2·log(rad_odd(p·q·N·|p-q|)) / log(2N)
  - ρ_proxy = 2·log(rad_odd(p·q·N)) / log(2N)
  - gap δ = ρ_true - ρ_proxy

Then study: mean, variance, extrema as functions of N,
and connect to Hardy-Littlewood prediction.
"""
import math
import json
import time
from collections import defaultdict

# ═══════════════════════════════════════════════════════════════
# PHASE 0: Prime sieve and utility functions
# ═══════════════════════════════════════════════════════════════
LIMIT = 20002

def sieve_primes(limit):
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime

print("Sieving primes...", end=" ", flush=True)
t0 = time.time()
IS_PRIME = sieve_primes(LIMIT)
PRIMES = [p for p in range(2, LIMIT + 1) if IS_PRIME[p]]
PRIME_SET = set(PRIMES)
print(f"done ({len(PRIMES)} primes up to {LIMIT}, {time.time()-t0:.2f}s)")

def rad_odd(n):
    """Product of distinct odd prime factors of n."""
    if n == 0:
        return 1
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    if n <= 1:
        return 1
    result = 1
    d = 3
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            result *= d
            while temp % d == 0:
                temp //= d
        d += 2
    if temp > 1:
        result *= temp
    return result

def odd_prime_set(n):
    """Set of odd prime factors."""
    if n == 0:
        return set()
    n = abs(n)
    factors = set()
    while n % 2 == 0:
        n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        factors.add(n)
    return factors

# ═══════════════════════════════════════════════════════════════
# PHASE 1: Enumerate all Goldbach pairs and compute ρ
# ═══════════════════════════════════════════════════════════════
print("\nPHASE 1: Enumerating Goldbach pairs for 2N = 6..10000")
print("=" * 60)

# For each even 2N, store: list of (p, q, rho_true, rho_proxy, delta)
results = {}
t0 = time.time()

for twoN in range(6, 10002, 2):
    N = twoN // 2
    log_twoN = math.log(twoN)
    rad_N = rad_odd(N)
    
    pairs = []
    for p in PRIMES:
        if p >= N:
            break
        q = twoN - p
        if q > p and q <= LIMIT and IS_PRIME[q]:
            # Both p, q are odd primes with p < q, p + q = 2N
            if p == 2:
                continue  # skip p=2 (even)
            
            diff = q - p  # always positive since q > p
            
            # ρ_proxy = 2·log(rad_odd(p)·rad_odd(q)·rad_odd(N)) / log(2N)
            # Since p, q are primes: rad_odd(p) = p, rad_odd(q) = q
            proxy_rad = p * q * rad_N
            rho_proxy = 2 * math.log(proxy_rad) / log_twoN
            
            # ρ_true = 2·log(rad_odd(p·q·N·(p-q))) / log(2N)
            # = 2·log(p · q · rad_odd(N) · rad_odd_new(diff)) / log(2N)
            known_primes = odd_prime_set(p) | odd_prime_set(q) | odd_prime_set(N)
            diff_primes = odd_prime_set(diff)
            new_primes = diff_primes - known_primes
            
            true_rad = proxy_rad
            for r in new_primes:
                true_rad *= r
            
            rho_true = 2 * math.log(true_rad) / log_twoN
            delta = rho_true - rho_proxy
            
            pairs.append((p, q, rho_true, rho_proxy, delta))
    
    if pairs:
        rhos_true = [x[2] for x in pairs]
        rhos_proxy = [x[3] for x in pairs]
        deltas = [x[4] for x in pairs]
        
        results[twoN] = {
            'N': N,
            'count': len(pairs),
            'rho_true_mean': sum(rhos_true) / len(rhos_true),
            'rho_true_std': (sum((x - sum(rhos_true)/len(rhos_true))**2 
                            for x in rhos_true) / len(rhos_true))**0.5,
            'rho_true_min': min(rhos_true),
            'rho_true_max': max(rhos_true),
            'rho_proxy_mean': sum(rhos_proxy) / len(rhos_proxy),
            'delta_mean': sum(deltas) / len(deltas),
            'delta_max': max(deltas),
        }

elapsed = time.time() - t0
print(f"Processed {len(results)} even numbers in {elapsed:.1f}s")
print(f"Total Goldbach pairs found: {sum(r['count'] for r in results.values())}")

# ═══════════════════════════════════════════════════════════════
# PHASE 2: Hardy-Littlewood comparison
# ═══════════════════════════════════════════════════════════════
print("\nPHASE 2: Hardy-Littlewood prediction")
print("=" * 60)

# Twin prime constant C₂ = ∏_{p>2} (1 - 1/(p-1)²)
C2 = 1.0
for p in PRIMES:
    if p > 2 and p < 1000:
        C2 *= (1 - 1.0 / (p - 1)**2)

print(f"C₂ (twin prime constant) ≈ {C2:.6f}")

def hardy_littlewood(twoN):
    """G(2N) ~ 2·C₂·∏_{p|N,p>2} (p-1)/(p-2) · N/ln(N)²"""
    N = twoN // 2
    if N < 3:
        return 0
    product = 1.0
    temp = N
    d = 3
    while d * d <= temp:
        if temp % d == 0:
            product *= (d - 1) / (d - 2)
            while temp % d == 0:
                temp //= d
        d += 2
    if temp > 2:
        product *= (temp - 1) / (temp - 2)
    
    return 2 * C2 * product * N / (math.log(N))**2

# Compare actual vs HL prediction for select values
print(f"\n{'2N':>6} {'G(2N)':>6} {'HL pred':>8} {'ratio':>7}")
print("-" * 35)
for twoN in [100, 200, 500, 1000, 2000, 5000, 10000]:
    if twoN in results:
        actual = results[twoN]['count']
        predicted = hardy_littlewood(twoN)
        ratio = actual / predicted if predicted > 0 else 0
        print(f"{twoN:>6} {actual:>6} {predicted:>8.1f} {ratio:>7.3f}")

# ═══════════════════════════════════════════════════════════════
# PHASE 3: Decomposition of ⟨ρ⟩ into static + dynamic
# ═══════════════════════════════════════════════════════════════
print("\nPHASE 3: Conductor ratio decomposition")
print("=" * 60)

# ρ_true = 2·[log(p) + log(q) + log(rad_N) + log(rad_new(|p-q|))] / log(2N)
#        = ρ_static(N) + ρ_boundary(p,q) + ρ_gap(p,q)
#
# where:
#   ρ_static  = 2·log(rad_odd(N)) / log(2N)         [depends only on N]
#   ρ_boundary = 2·[log(p) + log(q)] / log(2N)       [depends on pair]
#   ρ_gap     = 2·log(rad_odd_new(|p-q|)) / log(2N)  [noise term]

decomp_data = []
for twoN in sorted(results.keys()):
    if twoN < 20:
        continue
    N = twoN // 2
    log_twoN = math.log(twoN)
    
    rho_static = 2 * math.log(max(rad_odd(N), 1)) / log_twoN
    rho_mean = results[twoN]['rho_true_mean']
    delta_mean = results[twoN]['delta_mean']
    rho_proxy_mean = results[twoN]['rho_proxy_mean']
    count = results[twoN]['count']
    std = results[twoN]['rho_true_std']
    
    # ρ_boundary_mean = ρ_proxy_mean - ρ_static
    rho_boundary_mean = rho_proxy_mean - rho_static
    
    decomp_data.append({
        'twoN': twoN,
        'N': N,
        'count': count,
        'rho_mean': rho_mean,
        'rho_std': std,
        'rho_static': rho_static,
        'rho_boundary_mean': rho_boundary_mean,
        'delta_mean': delta_mean,
        'log2N': math.log(twoN),
    })

# ═══════════════════════════════════════════════════════════════
# PHASE 4: Key statistics and asymptotic fits
# ═══════════════════════════════════════════════════════════════
print("\nPHASE 4: Asymptotic analysis")
print("=" * 60)

# Extract data for fitting (use 2N >= 100 for asymptotic regime)
fit_data = [d for d in decomp_data if d['twoN'] >= 100]

log2N_vals = [d['log2N'] for d in fit_data]
rho_means = [d['rho_mean'] for d in fit_data]
rho_stds = [d['rho_std'] for d in fit_data]
rho_statics = [d['rho_static'] for d in fit_data]
delta_means = [d['delta_mean'] for d in fit_data]
counts = [d['count'] for d in fit_data]

# Linear regression: ρ_mean = a·log(2N) + b
n = len(fit_data)
sx = sum(log2N_vals)
sy = sum(rho_means)
sxx = sum(x**2 for x in log2N_vals)
sxy = sum(x*y for x, y in zip(log2N_vals, rho_means))
a_rho = (n * sxy - sx * sy) / (n * sxx - sx**2)
b_rho = (sy - a_rho * sx) / n

# R² for the fit
ss_res = sum((y - (a_rho * x + b_rho))**2 for x, y in zip(log2N_vals, rho_means))
ss_tot = sum((y - sy/n)**2 for y in rho_means)
r2_rho = 1 - ss_res / ss_tot

print(f"⟨ρ_true⟩ ≈ {a_rho:.4f} · log(2N) + {b_rho:.4f}")
print(f"R² = {r2_rho:.6f}")
print()

# Fit for ⟨δ⟩ (gap term)
sy_d = sum(delta_means)
sxy_d = sum(x*y for x, y in zip(log2N_vals, delta_means))
a_delta = (n * sxy_d - sx * sy_d) / (n * sxx - sx**2)
b_delta = (sy_d - a_delta * sx) / n

print(f"⟨δ⟩ ≈ {a_delta:.4f} · log(2N) + {b_delta:.4f}")

# Fit for σ(ρ)
sy_s = sum(rho_stds)
sxy_s = sum(x*y for x, y in zip(log2N_vals, rho_stds))
a_std = (n * sxy_s - sx * sy_s) / (n * sxx - sx**2)
b_std = (sy_s - a_std * sx) / n

print(f"σ(ρ) ≈ {a_std:.4f} · log(2N) + {b_std:.4f}")

# Key insight: what fraction of ρ comes from static vs dynamic?
print(f"\nDecomposition at 2N = 10000:")
d10k = [d for d in decomp_data if d['twoN'] == 10000][0]
total = d10k['rho_mean']
static_frac = d10k['rho_static'] / total * 100
boundary_frac = d10k['rho_boundary_mean'] / total * 100
gap_frac = d10k['delta_mean'] / total * 100
print(f"  ⟨ρ_true⟩    = {total:.4f}")
print(f"  ρ_static    = {d10k['rho_static']:.4f}  ({static_frac:.1f}%)")
print(f"  ⟨ρ_boundary⟩ = {d10k['rho_boundary_mean']:.4f}  ({boundary_frac:.1f}%)")
print(f"  ⟨δ⟩          = {d10k['delta_mean']:.4f}  ({gap_frac:.1f}%)")
print(f"  σ(ρ)        = {d10k['rho_std']:.4f}")
print(f"  G(10000)    = {d10k['count']}")

# ═══════════════════════════════════════════════════════════════
# PHASE 5: The conductor density — a new quantity
# ═══════════════════════════════════════════════════════════════
print("\nPHASE 5: Conductor density Σ(N)")
print("=" * 60)

# Define: Σ(N) = Σ_{p+q=2N} ρ_true(p,q) / G(2N)²
# This measures the "total conductor weight per pair"
# If it has a clean asymptotic, it connects ρ to G(2N)

print(f"\n{'2N':>6} {'G(2N)':>6} {'⟨ρ⟩':>8} {'σ(ρ)':>7} "
      f"{'⟨ρ⟩/log(2N)':>11} {'σ/⟨ρ⟩':>7}")
print("-" * 55)
for twoN in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
    if twoN in results:
        r = results[twoN]
        ratio = r['rho_true_mean'] / math.log(twoN)
        cv = r['rho_true_std'] / r['rho_true_mean'] if r['rho_true_mean'] > 0 else 0
        print(f"{twoN:>6} {r['count']:>6} {r['rho_true_mean']:>8.4f} "
              f"{r['rho_true_std']:>7.4f} {ratio:>11.4f} {cv:>7.4f}")

# ═══════════════════════════════════════════════════════════════
# SAVE DATA for plotting
# ═══════════════════════════════════════════════════════════════
print("\nSaving data...")

# Save summary for every 2N
with open('/home/claude/paper14/data/conductor_census.json', 'w') as f:
    # Convert to list for JSON
    out = []
    for twoN in sorted(results.keys()):
        r = results[twoN]
        r['twoN'] = twoN
        out.append(r)
    json.dump(out, f)

# Save decomposition data
with open('/home/claude/paper14/data/decomposition.json', 'w') as f:
    json.dump(decomp_data, f)

# Save fit parameters
fit_params = {
    'rho_mean': {'a': a_rho, 'b': b_rho, 'R2': r2_rho},
    'delta_mean': {'a': a_delta, 'b': b_delta},
    'rho_std': {'a': a_std, 'b': b_std},
    'C2': C2,
}
with open('/home/claude/paper14/data/fit_params.json', 'w') as f:
    json.dump(fit_params, f, indent=2)

print("All data saved.")
print(f"\nTotal Goldbach pairs processed: {sum(r['count'] for r in results.values())}")
