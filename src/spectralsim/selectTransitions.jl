struct TransitionsInc
    a::Vector{Float32}   # amplitudes
    γ::Vector{Float32}   # linewidths
    ωᵣ::Vector{Float32}  # transition frequencies
    branch::Vector{String}  # branch of transitions
    v::Vector{Int32}      # ground state v   
    J::Vector{Int32}      # ground state J
end

function selectTransitions(transitions,γ,amplitudes,ROI)
    #####
    ## THIS IS FAR FROM OPTIMAL AS IT CONVERTS BACK AND FORTH BETWEEN FLOAT64 AND INT
    #####
#     for specimen = s.species
#         % transitions outside range have been blanked with NaN
#         % serialize all transitions with ground state J and V. Also store a
#         % flag for the transition (-2 for O, 0 for Q, +2 for S)
#         % also, if dual pump is enabled by a non-zero frequency offset, store
#         % the pumping region (1 or 2). this is important for the convolution to
#         % be performed later.
         J_extended = repeat(J,length(V[1:end-1]))
         V_extended = vec(reshape(repeat(V[1:end-1]',length(J)),(:,1)))

         O = fill!(similar(J_extended),-2)
         Q = fill!(similar(J_extended),0)
         S = fill!(similar(J_extended),2)
         PR = fill!(Vector{Float64}(undef,length(J_extended)),NaN) # Is in Pumping region? preallocated as NaN

         transitions = [vec(reshape(transitions.O,(:,1))) V_extended J_extended O PR;
                        vec(reshape(transitions.Q,(:,1))) V_extended J_extended Q PR;
                        vec(reshape(transitions.S,(:,1))) V_extended J_extended S PR];

    
#         % decide if a transition falls into pumping region. assign NaN
#         % if totally outside the plotting region. these are excluded from
#         % further computations
         pumpregion1 = @. (transitions[:,1]>ROI[1]) && (transitions[:,1]<ROI[2])
         transitions[pumpregion1,5].=1;

#         % remove all transitions that are outside pumping region (+ some margin)
        transitions = transitions[vec(.!any(isnan.(transitions),dims=2)),:]

        # throw an error if it is outside the region
        if isempty(transitions)
            error("No transitions in selected ROI, aborting.")
        end

        # extract relevant data
        ωᵣ = transitions[:,1]

        # extract the relevant J-dependend linewidths
        γᵢ = γ[Int.(transitions[:,3].+1)]

        # extract the amplitudes based on branch, v and J
        a = similar(γᵢ)
        v = Vector{Int32}(undef,length(a))
        Js = Vector{Int32}(undef,length(a))
        branch = Vector{String}(undef,length(a))
        for i = 1:size(transitions,1)
            Js[i],v[i] = Int(transitions[i,3]), Int(transitions[i,2])
            if transitions[i,4] == -2   # O branch
                a[i] = amplitudes.O[Js[i]+1,v[i]+1]
                branch[i] = "O"
            elseif transitions[i,4] == 0   # Q branch
                a[i] = amplitudes.Q[Js[i]+1,v[i]+1]
                branch[i] = "Q"
            elseif transitions[i,4] == 2   # S branch
                a[i] = amplitudes.S[Js[i]+1,v[i]+1]
                branch[i] = "S"
            end
            
        end

        # only include transitions with > 0.001 of max amplitude
        idx=a.>0.001*maximum(a)

        # store in output struct and convert to Float32 at this point
        # would probably be better to calculate it with single precision in the first place
        # alright for now.
        transitions_inc = TransitionsInc(Float32.(a[idx]),Float32.(γᵢ[idx]),Float32.(ωᵣ[idx]),branch[idx],v[idx],Js[idx])

return transitions_inc
end