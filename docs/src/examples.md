# Examples

## Simulation of Single-Pump Theoretical Susceptibility Spectra
```@example
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
savefig(p, "figure1.html") # ignore this line, this is just for the documentation build process
```

```@raw html
<iframe src="../figure1.html" style="height:500px;width:100%;"></iframe>
```