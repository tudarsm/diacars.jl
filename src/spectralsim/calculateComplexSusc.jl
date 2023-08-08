function isolatedLines!(χR::Vector{Float32}, χI::Vector{Float32}, a::Vector{Float32}, ωᵣ::Vector{Float32}, ω::Vector{Float32}, γᵢ::Vector{Float32})
    @avxt for i = eachindex(γᵢ)
        for j = eachindex(ω)
            χI[j] += a[i] * γᵢ[i] / (γᵢ[i]^2 + (ωᵣ[i] - ω[j])^2)
        end
    end
    @avxt for i = eachindex(γᵢ)
        for j = eachindex(ω)
            χR[j] += a[i] * (ωᵣ[i] - ω[j]) / (γᵢ[i]^2 + (ωᵣ[i] - ω[j])^2)
        end
    end
    nothing
end

function voigtLines!(χR::Vector{Float32}, χI::Vector{Float32}, a::Vector{Float32}, ωᵣ::Vector{Float32}, ω::Vector{Float32}, γᵢ::Vector{Float32},MolecularParameters,ωres,T)
    # need two more scratch arrays
    χRc = similar(χR)
    χIc = similar(χI)
    fill!(χRc,0)
    fill!(χIc,0)

    # iterate through transitions
    for i = eachindex(γᵢ)
        # construct a Gaussian kernel with appropriate size and width for the current transition
        γd = convert(Float32,4.3014e-7*ωᵣ[i]*sqrt(T/MolecularParameters.CMASS))  # gaussian width
        Γd = gaussiankernelT(ωres,γd)                           # gaussian kernel

        # generate the current isolated line transitions
        @avxt for j = eachindex(ω)
            χRc[j] = a[i] * (ωᵣ[i] - ω[j]) / (γᵢ[i]^2 + (ωᵣ[i] - ω[j])^2)
        end
        @avxt for j = eachindex(ω)
            χIc[j] = a[i] * γᵢ[i] / (γᵢ[i]^2 + (ωᵣ[i] - ω[j])^2)
        end

        # now convolve the results with the Doppler width and sum up
      
        χR .+= χRc⋆Γd
        χI .+= χIc⋆Γd

    end

    nothing
end

function rotationalDiffusion!(χR::Vector{Float32}, χI::Vector{Float32}, a::Vector{Float32}, ωᵣ::Vector{Float32}, ω::Vector{Float32}, γᵢ::Vector{Float32},fBrot,γⱼ,P,TT,transitions_inc,ωres,T,MolecularParameters,usedoppler)

    # get the denominator for rotational diffusion
    denominator = getDenominatorForRotationalDiffusion(fBrot,P,γⱼ,TT,ω)
   
    # iterate through transitions
    # need a complex scratch array for the current transition
    χc = Vector{ComplexF32}(undef,length(χR))
    fill!(χc,0)

    for i = eachindex(γᵢ)
        # calculate isolated lines for this transition
        @. χc = a[i]/(ωᵣ[i]-ω-im*γᵢ[i])

        # if the underlying line shape should be Voigt, do convolution with Doppler width here
        if usedoppler
            γd = convert(Float32,4.3014e-7*ωᵣ[i]*sqrt(T/MolecularParameters.CMASS))  # gaussian width
            Γd = gaussiankernelT(ωres,γd)                           # gaussian kernel
            # reinterpret the complex array as two independent arrays
            y = reinterpret(real(eltype(χc)), χc)
            y = reshape(y, 2, size(χc)...)
            y[1,:] .= y[1,:]⋆Γd
            y[2,:] .= y[2,:]⋆Γd
        end

        # if it is a q-branch transition, use rotational diffusion. o and s are assumed to be isolated
        if transitions_inc.branch[i] == "Q"
            @. χc = χc * denominator[transitions_inc.v[i]+1,:]
        end

        # add to χR and χI
        @. χR += real(χc)
        @. χI += imag(χc)
    end


end


function getDenominatorForRotationalDiffusion(fBrot,P,γⱼ,TT,ω)

        # PORTED FROM MATLAB
    # We need the denominator from equation (4) in Hall and Greenhalgh.
    # the numerator and the term before the fraction is just the
    # ordinary lorentzian equation for chi(3). basically, the
    # transitions are weighted with the factor 1/(1+i<[f/tau]/Dv>)
    f = fBrot[:,1:end-1]                    # this is f in equation (4) 
    τⱼ = @. 1/(P*γⱼ);                       # this is τⱼ in equation(4)
    ωvvp1 = diff(TT,dims=2)                 # this is omega_v,v+1(omega_r) (i.e. all line positions of the q branches for each v)
    
    # this is not the fastest way probably, but it works.
    ωs = copy(ω)
    ωs = reshape(ωs,(1,1,length(ω)))        # add singleton dimensions
    ωvvp1r = repeat(ωvvp1,outer=[1,1,length(ω)])
    Dv = @. ωs-ωvvp1r - im/τⱼ
    denominator = dropdims(1 .+ im.*sum(f./τⱼ./Dv;dims=1);dims=1) # expectation value in J direction
    @. denominator = denominator/denominator^2;       # this is still somewhat unclear
    
    return denominator
end

function calculateComplexSusc(transitions_inc,Model,JMax,γⱼ,P,ROI,T,ωres,fBrot,MolecularParameters,TT)
    # construct wavenumberarray based on linewidth
    
    #ωres = round(0.1*minimum(γ[1:JMax])*P,sigdigits=1)
    
    ωres = convert(Float32,ωres)    # make it a Float32
    ω = convert(Vector{Float32}, collect(ROI[1]:ωres:ROI[2]))
    
    # preallocate χ with same length as ω
    χR = Vector{Float32}(undef, length(ω))
    χI = Vector{Float32}(undef, length(ω))
    # ...and fill it with 0
    fill!(χR,0)
    fill!(χI,0)

    N = Float32(c.NA * 273.15 / T * P / c.molvol)    # number density

    a = N.*transitions_inc.a
    ωᵣ = transitions_inc.ωᵣ
    γ = transitions_inc.γ/2*P     # divide by two here instead of in Lorentzian and don't forget to mulitply the pressure

    if Model == "I"
        isolatedLines!(χR,χI,a,ωᵣ,ω,γ)
    elseif Model == "R"
        rotationalDiffusion!(χR,χI,a,ωᵣ,ω,γ,fBrot,γⱼ,P,TT,transitions_inc,ωres,T,MolecularParameters,false)
    elseif Model == "V"
        voigtLines!(χR,χI,a,ωᵣ,ω,γ,MolecularParameters,ωres,T)
    elseif Model == "VR"
        rotationalDiffusion!(χR,χI,a,ωᵣ,ω,γ,fBrot,γⱼ,P,TT,transitions_inc,ωres,T,MolecularParameters,true)
    else
        error("Model has to be either (I)solaten lines, (R)otational Diffusion, (V)oigt, or Rotational Diffusion with Voigt (VR).")
    end

    return χR,χI,ω,ωres

end