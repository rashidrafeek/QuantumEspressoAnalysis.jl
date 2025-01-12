# QuantumEspressoAnalysis

[![Build Status](https://github.com/rashidrafeek/QuantumEspressoAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/rashidrafeek/QuantumEspressoAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)

# QuantumEspressoAnalysis

A Julia module to simplify the analysis of Quantum Espresso (QE) output files.

## Example usage

```julia
import QuantumEspressoAnalysis as QE

energy = QE.extract_final_energy("output_scf_file.out")
println("Final energy: $energy")
```
