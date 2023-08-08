using LinearAlgebra
using Plots
#using BenchmarkTools
using LoopVectorization
using FourierTools
using DSP

#### ROADMAP / TO DOS
# create git repo
# clean up functions, add references and meaningful comments
# optimize functions
# allow multi-species / single-pump and dual-pump
# library generation
# gui?
# documentation
# examples
# benchmarks / tests
# multi threading
# library based spectra generation
# mixed-integer fitting


include("getMolecularParameters.jl")
include("calculateTermEnergies.jl")
include("calculateLinePositions.jl")
include("calculateLineWidths.jl")
include("calculateBoltzmannFraction.jl")
include("constants.jl")
include("chiamplitudes.jl")
include("selectTransitions.jl")
include("calculateComplexSusc.jl")

#function doAll(T)
    T=2000
    Species = "N2"
    X = 0.5
    P = 1
    CHINR_Buffer = 10
    Model='I'
    ROI = [2200 2400]

    N2Param = getMolecularParameters(Species)
    TT,FJ,GJ,G0V = calculateTermEnergies(N2Param)
    transitions = calculateLinePositions(TT)
    γ = calculateLineWidths(T,FJ,Species)
    Δ,JMax,fBrot = calculateBoltzmannFraction(T,FJ,GJ,G0V)
    amplitudes = chiamplitudes(N2Param,Δ)
    transitions_inc = selectTransitions(transitions,γ,amplitudes,ROI)
    @time χR,χI,ω,ωres = calculateComplexSusc(transitions_inc,Model,JMax,γ,P,ROI,T)
    @time χR,χI,ω,ωres = calculateComplexSusc(transitions_inc,Model,JMax,γ,P,ROI,T)
    display(plot(ω,χR))
    # println("================TOTAL================")
    #return ω,I
    
    # for the convolution, check out FourierTools.jl
    # https://bionanoimaging.github.io/FourierTools.jl/dev/convolutions/#FourierTools.plan_conv


#end

#ω,I=doAll(1600);
#display(plot(ω,I))
# doAll(2000)
# doAll(2000)
#println("TOTAL")
#@time for T = 300:50:2200
#@time for T = 300:50:2200
#    ω,I=doAll(T)
#    display(plot(ω,I))
#end
# println("================TIMING RESULTS================")
# @time for i=1:100
#  doAll(2000)
# end
# println("================TIMING RESULTS================")
# plot(χ)