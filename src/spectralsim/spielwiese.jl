using DiaCARS
using Plots

Species = "CO"
T=300

#
P=1
Model ="I"
ROI=[2100 2435]
ωres = 1e-3

@time χR,χI,ω = simulateTheoreticalSusceptibility(T,Species,P,Model,ROI,ωres)

χres = χI.^2 .+ χR.^2
@show maximum(sqrt.(χres))
plot(ω,sqrt.(χres))
# ylims!(0,1000)
