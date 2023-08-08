module DiaCARS

using LinearAlgebra
using LoopVectorization
using FourierTools
using DSP
using FFTW

# include files
include("spectralsim/getMolecularParameters.jl")
include("spectralsim/calculateTermEnergies.jl")
include("spectralsim/calculateLinePositions.jl")
include("spectralsim/calculateLineWidths.jl")
include("spectralsim/calculateBoltzmannFraction.jl")
include("spectralsim/constants.jl")
include("spectralsim/chiamplitudes.jl")
include("spectralsim/selectTransitions.jl")
include("spectralsim/calculateComplexSusc.jl")

export simulateTheoreticalSusceptibility
export simulateConvolvedSpectrum
export calculateBoltzmannFraction

"""
    simulateTheoreticalSusceptibility(T,Species,P,Model,ROI,FrequencyOffset)

Simulates the theoretical suscepbtility of a given Species based on temperature T in K, pressure P in bar on a wavenumberarray ω based on the region of interest ROI 
The calculated amplitudes are valid only for single-pump CARS. For Dual-Pump Spectra, the results have to be divided by two.

- T: Temperature in K
- Species: "N2","O2" or "CO"
- P: Pressure in Bar
- Model: "I" for isolated lines, "V" for Voigt profile or "R" for rotational diffusion
- ROI: [ωmin ωmax], range for which to calculate the suscepbtility. Expect edge effects if using Voigt
- ωres: Resolution of original wavenumber array. Should be sufficiently small to resolve transitions.
"""
function simulateTheoreticalSusceptibility(T,Species,P,Model,ROI,ωres)

    MolecularParameters = getMolecularParameters(Species)
    TT,FJ,GJ,G0V = calculateTermEnergies(MolecularParameters)
    transitions = calculateLinePositions(TT)
    γ = calculateLineWidths(T,FJ,Species)
    Δ,JMax,fBrot = calculateBoltzmannFraction(T,FJ,GJ,G0V)
    amplitudes = chiamplitudes(MolecularParameters,Δ)
    transitions_inc = selectTransitions(transitions,γ,amplitudes,ROI)
    χR,χI,ω = calculateComplexSusc(transitions_inc,Model,JMax,γ,P,ROI,T,ωres,fBrot,MolecularParameters,TT)
    return χR,χI,ω,Species,MolecularParameters

end


function simulateConvolvedSpectrum(χR,χI,γ,method,ω)
    I = convolveTheoSusc(χR,χI,1,"single",ω,ωres)
    return I
end


⋆(a,b) = dconv_edge(a,b)

function dconv_edge(signal::Vector{T},kern::Vector{T}) where {T}
    M = length(kern)
    N = length(signal)
    m = ceil(Int,M/2) # this has to be ceil in case the length of the kernel is odd.
    out = zeros(T,N)
    @avxt for j in M:N
        tmp = zero(T)
        for i in 1:M
            tmp += signal[j-i+1] * kern[i]
        end
        out[j-m+1]=tmp
    end
    return out
end



function gaussiankernelT(ωres::T,wid::T) where {T<:Real}
    # get a grid 4 times the width
    # at this width, g[1]/maximum(g) ≈ 1.5e-5 which should be close enough
    x = collect(T,0:ωres:2*wid)
    # force it to be symmetric. not necessarily true if -2:ωres:2
    x = vcat(reverse(x[2:end]),x)
    #x = collect(T,-2:ωres:2)
    # get gaussian
    c = convert(T,0.6006)
    g = @. exp(-((x)./(c.*wid)) .^2);
    # normalize
    g = g./sum(g);
    return g
end

end
