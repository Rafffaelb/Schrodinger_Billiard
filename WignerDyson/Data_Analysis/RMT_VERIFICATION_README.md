# RMT Verification Tools

This directory contains tools to verify that your Hamiltonian construction correctly implements Random Matrix Theory (RMT) for the three Wigner-Dyson ensembles.

## Quick Start

### Option 1: Jupyter Notebook (Recommended for detailed analysis)

```bash
jupyter notebook RMT_Verification.ipynb
```

The notebook includes:
- Interactive verification of all three ensembles (GOE, GUE, GSE)
- Level spacing distribution comparisons
- Variance analysis of matrix elements
- Beautiful visualizations
- Statistical tests (Kolmogorov-Smirnov)

### Option 2: Python Script (Quick command-line verification)

```bash
# Verify all ensembles with default parameters
python rmt_verify.py

# Verify only GOE
python rmt_verify.py --ensemble GOE

# Custom parameters
python rmt_verify.py --ensemble GUE --matrix-size 200 --samples 1000 --variance 1.0

# See all options
python rmt_verify.py --help
```

## What Gets Verified?

### 1. Level Spacing Distribution
The **gold standard** test for RMT. Compares your Hamiltonian's eigenvalue spacings against the theoretical Wigner surmise:

- **GOE**: P(s) = (π/2)s exp(-πs²/4)
- **GUE**: P(s) = (32/π²)s² exp(-4s²/π)
- **GSE**: P(s) = (2¹⁸/3⁶π³)s⁴ exp(-64s²/9π)

### 2. Variance of Matrix Elements
Checks that diagonal and off-diagonal elements have the correct variance according to RMT normalization.

### 3. Level Repulsion
Verifies that small spacings are suppressed (quantum chaos signature).

## Interpreting Results

### ✅ Success Indicators

1. **Histogram matches Wigner surmise** - Your Hamiltonian construction is correct!
2. **KS statistic < 0.05** - Statistical test confirms the match
3. **Mean spacing ≈ 1** - Proper spectrum unfolding
4. **Few spacings < 0.1** - Level repulsion is present

### ⚠️ Warning Signs

1. **Histogram doesn't match** - Check variance normalization
2. **KS statistic > 0.1** - Significant deviation from theory
3. **Many spacings near 0** - Missing level repulsion (wrong ensemble?)
4. **Variance values off** - Scaling factors in Hamiltonian construction need adjustment

## Understanding the Ensembles

| Ensemble | Symmetry | Matrix Type | Physical System |
|----------|----------|-------------|----------------|
| **GOE** | Real symmetric | H = H^T | Time-reversal symmetric, spinless |
| **GUE** | Complex Hermitian | H = H† | Time-reversal broken |
| **GSE** | Quaternion self-dual | Special structure | Time-reversal symmetric with spin-orbit |

## Output Files

The script generates PNG files with visualizations:
- `GOE_spacing_distribution.png`
- `GUE_spacing_distribution.png`
- `GSE_spacing_distribution.png`

Each shows:
- **Left panel**: Histogram vs Wigner surmise
- **Right panel**: Cumulative distribution comparison

## Dependencies

```bash
pip install numpy scipy matplotlib seaborn jupyter
```

## Troubleshooting

### Problem: Distributions don't match

**Solution**: Try different variance values:
```bash
python rmt_verify.py --variance 0.5
python rmt_verify.py --variance 2.0
```

The level spacing distribution is **universal** and should match regardless of variance, but incorrect variance might cause numerical issues.

### Problem: "Variance values are wrong"

This indicates the diagonal/off-diagonal elements don't have the expected variance. You may need to adjust the scaling factors in your C++ Hamiltonian construction:

```cpp
// In Orthogonal.cpp, Unitary.cpp, Symplectic.cpp
// Adjust the sqrt(V/2.0) and sqrt(V/1.0) factors
```

Refer to the detailed analysis in `hamiltonian_analysis.md` in the brain artifacts.

## References

1. Mehta, M.L. "Random Matrices" (2004)
2. Haake, F. "Quantum Signatures of Chaos" (2010)
3. Beenakker, C.W.J. "Random-matrix theory of quantum transport" Rev. Mod. Phys. 69, 731 (1997)

## Next Steps After Verification

If your verification passes ✅:
1. Your Hamiltonian construction is correct!
2. You can confidently run your quantum chaos simulations
3. Move on to analyzing conductance, entanglement, etc.

If verification fails ❌:
1. Check the variance analysis output
2. Adjust scaling factors in your C++ code
3. Re-run verification
4. Consult `hamiltonian_analysis.md` for detailed guidance
