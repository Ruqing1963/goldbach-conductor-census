# Conductor Census of 425,082 Goldbach–Frey Curves

**Paper #14 in the Titan Project series**

Large-scale computation of conductor ratios for all Goldbach pairs with p + q ≤ 10,000, establishing asymptotic decomposition and bandwidth stability.

## Key Findings

1. **Three-component decomposition**: ρ = ρ_static (6.5%) + ρ_boundary (66.0%) + δ (27.5%)
2. **Bandwidth stability**: σ(ρ) ≈ 0.33 ± 0.05, constant across 2N ∈ [100, 10000]
3. **Signal-to-noise improvement**: CV = σ/⟨ρ⟩ → 0 as N → ∞
4. **Hardy–Littlewood consistency**: G(2N) tracks HL prediction with expected fluctuations
5. **Conjecture**: ⟨ρ⟩ = 4 + 2·log(rad_odd(N))/log(2N) + O(log log N / log N)

## Census Statistics

| 2N | G(2N) | ⟨ρ⟩ | σ(ρ) | CV |
|------|-------|-------|------|------|
| 100 | 6 | 4.992 | 0.349 | 7.0% |
| 1000 | 28 | 5.160 | 0.369 | 7.2% |
| 10000 | 127 | 5.377 | 0.323 | 6.0% |

Total: 425,082 Goldbach pairs across 4,997 even numbers.

## Repository Structure

```
├── paper/          Conductor_Census.pdf, .tex
├── figures/        3 figure pairs (PDF + PNG)
├── data/           Census JSON data files
├── scripts/        Computation and figure generation
├── README.md, LICENSE
```

## Series Context

| # | Paper | DOI |
|---|-------|-----|
| 12 | True Conductor Validation | 10.5281/zenodo.18749731 |
| 13 | Universal Tame Semistability | 10.5281/zenodo.18751169 |
| **14** | **Conductor Census** | *this paper* |

## License

MIT
