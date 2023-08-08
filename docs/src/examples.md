# Examples

## Simulation of Single-Pump Theoretical Susceptibility Spectra
```julia
using DiaCARS
using Plots

T=2000
Species = "N2"
ROI = [2200 2400]
FrequencyOffset = 0
Model = "I"
P = 1

χR,χI,ω = simulateTheoreticalSusceptibility(T,Species,P,Model,ROI,FrequencyOffset)

# Plot the result
plotlyjs()
p=plot(ω,[χR,χI],title="Theoretical Susceptibility", label=["Real" "Imaginary"])
xlabel!("Raman shift in cm-1")
ylabel!("Amplitude")
xlims!((ROI[1],ROI[2]))
```

