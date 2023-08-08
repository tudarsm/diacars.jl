struct Δs
    O
    Q
    S
end

function calculateBoltzmannFraction(T,FJ,GJ,G0V)
    # % calculate the differential boltzmann fraction used for the amplitudes of
    # % chi(3)
    # % max. vibrational state and vibrational partition function qv:
    # % a shortcut to hc/kT for easier readability.
    hc_KT = c.hc_k./T; # this is hc/kT
    

    #     % populations
    #     % the following looks a little complicated due to the vectorization
    #     % it just is (2J+1)*GJ*exp(-FJ) for every temperature in the range
    #     % (accounting for the different degeneracies for even and odd Js)
    #     % to do this vectorized, FJ is reshaped from NJxNV to NJ*NVx1 and the
    #     % result is reshaped back. See equation (III,161)
    #     % This is the numerator in the Boltzmann equation for all rotational
    #     % energy levels
    RS = @. (2 .* J.+1) .* GJ' .* exp.(-FJ*hc_KT)
    #     % This is the population distribution, accounting for the population of
    #     % the vibrational states
    #     VS(1,:,:) = exp(-s.mol.(specimen{:}).G0V'.*hc_KT); % add a singleton dimension for later repmat
    RSVS = @. RS * exp(-G0V'.*hc_KT); 
    #     % this is now the combination of vibrational and rotational states
    #     RSVS = repmat(VS,length(s.J),1,1).*RS;
    #     fB = RSVS./sum(sum(RSVS,1),2); % get the boltzmann distribution by normalizing to the sum of all states (i.e. the partition function)
    fB = RSVS./sum(RSVS[:])
    Q = -diff(fB,dims=2)
    O = similar(Q).*NaN
    S = similar(O)
    O[3:end,:] = -(fB[1:end-2,2:end]-fB[3:end,1:end-1,:]);
    S[1:end-2,:] = -(fB[3:end,2:end]-fB[1:end-2,1:end-1,:]);
    # fill!(S,0)
    # fill!(O,0)
    Δ = Δs(O,Q,S)
    
    #     % now, for every transition, get the differential population
    #     % distribution
    #     % Note: diff uses does i-(i+1), we want the opposite, so invert it.
    #     delta.Q = -diff(fB,1,2);
    #     delta.O = NaN(size(delta.Q));
    #     delta.O(3:end,:,:) = -(fB(1:end-2,2:end,:)-fB(3:end,1:end-1,:));
    #     delta.S = NaN(size(delta.Q));
    #     delta.S(1:end-2,:,:) = -(fB(3:end,2:end,:)-fB(1:end-2,1:end-1,:));
    
    #     % Extract the highest J level with significant population difference for the later
    #     % estimation of the desired grid size. This is useful because higher
    #     % J's tend to have lower linewidth. But if they are not populated, then
    #     % it doesn't make sense to resolve them.
    #     % Threshold is 0.001*max, i.e. 0.1% contribution. This does not mean that
    #     % the lines with lower values are disregarded, only that their
    #     % resolution is lower. limit to max s.J+1 (due to matlab starting arrays at index 1).
    #     s.JMax=min([round(1.5*find((delta.Q(:,1)./max(delta.Q(:,1)))>0.001,1,'last')) max(s.J)+1]);
    JMax = Int(minimum([round(1.5.*findlast(Q[:,1]./maximum(Q[:,1]).>0.001)) maximum(J)]))
    #     % rotational partition function, used later for rotational diffusion
    #     QR = s.mol.(specimen{:}).GJ.*(2*s.J+1)*exp(-s.mol.(specimen{:}).FJ(:,1)*hc_KT); % rotational partition function for every temperature
    QR =  (vec(GJ).*(2J.+1))'*exp.(-FJ[:,1]*hc_KT)
    fBrot = RS./QR

    #     % store for later use
    #     s.delta.(specimen{:}).Q = delta.Q;
    #     s.delta.(specimen{:}).O = delta.O;
    #     s.delta.(specimen{:}).S = delta.S;
    #     s.fB_rot = RS/QR; % for rotational diffusion. this is a normalized population distribution without accounting for the population of the vibrational band
    return Δ,JMax,fBrot,fB
    end