function calculateLineWidths(T,FJ,species)
    γ = lineWidthsFromMEG(T,FJ,species)
    return γ
end


function lineWidthsFromMEG(T,FJ,species)
    # % [1] From Koszykowski, Rahn, Palmer, Coltrin: Theoretical and
    #         % Experimental Studies of High-Resolution Inverse Raman Spectra
    #         % of N2 at 1-10 atm, Phys Chem 91,1,1987
    #         % AND
    #         % [2] Rahn, Palmer, Koszykowski, Greenhalgh: Comparison of rotationally
    #         % inelastic collision models for q-branch raman spectra of N2,
    #         % Chemical Physical Letters 133,6,1987. Primary source probably
    #         % L. A. Rahn and R. E. Palmer, "Studies of nitrogen self-broadening at high temperature with inverse Raman spectroscopy," J. Opt. Soc. Am. B 3, 1164-1169 (1986)
    #         % AND
    #         % [3] Sitz, Greg O., and R. L. Farrow. "Pump probe measurements of state to state rotational energy transfer rates in N2 (v= 1)." The Journal of chemical physics 93.11 (1990): 7883-7893.

    #         % some handles for better readability
    #         % data from [2]
    if species == "N2"
        α = 0.0231
        tcorr = sqrt(295 / T).*(1-exp(-0.1487))./(1-exp(-0.1487*T/295)) # temperature correction, [2] eq (1)
        β = 1.67
        δ = 1.21  
        a = 1.5
    elseif species == "O2"
        α = 0.0167
        tcorr = (295 / T).^(1.31); # with m = inf -> exp terms -> 1
        β = 1.45
        δ = 1.32
        a = 1.5
    elseif species == "CO"          # https://doi.org/10.1063/1.455764. Coefficients are slightly different than in CARSFT.
        α = 0.01334
        tcorr = sqrt(295 / T).*(1-exp(-0.19))./(1-exp(-0.19*T/295)) # temperature correction, [2] eq (1)
        β = 1.452
        δ = 1.246
        a = 2 
    end

    Eᵢ = FJ[:,1].* c.h .* c.c
    # make sure it is a length(J)x1 matrix
    Eᵢ = reshape(Eᵢ,length(Eᵢ),1)
    Eⱼ = FJ[:,1]' .* c.h .* c.c
    dEᵢⱼ = Eⱼ .- Eᵢ # this is \δ E_{ij}
    kT = c.k*T; # this is kT



    #         % this is without the pressure, this will be multiplied later
    #         % on so it is not included in the linewidth calculation yet
    #         % the transponings are necessary to account for the switching
    #         % between ij and ji.
        γⱼᵢ = tril( transpose((α*tcorr*((1 .+ a .* Eᵢ ./ (kT.*δ)) ./ (1 .+ a.*Eᵢ./kT)).^2 .* exp.(-β.*dEᵢⱼ./kT))),-1)
    #         % this is microscopic reversibility [1], eq (4.2):
    γᵢⱼ = triu(
            transpose(
                (
                    transpose(2 .*J .+ 1)./(2 .* J .+ 1)
                    .*
                    γⱼᵢ
                    .*
                    transpose(exp.(dEᵢⱼ ./kT))
                )
            )
        ,1);
            #         % combine the upper and lower triangular matrices to the entire
    #         % S matrix
             S = γⱼᵢ + γᵢⱼ;
             
    #         % account for selection rules DJ = +-2
    #         % this basically means, that Ji and Jj are Eᵢther both even or
    #         % uneven. this introduces a checkerboard like pattern. this is
    #         % multiplied onto the calculated S matrix to make the
    #         % forbidden entries zero.
    # DEACTIVATED FOR CO
             checkerboard = ones(length(J),length(J))
             if species != "CO"
                checkerboard[2:2:(J[end]+1)^2].=0;
             end
             S = S .* checkerboard;

    #         % this is equation [2] (4.3). take the upper triangular matrix of
    #         % γᵢⱼ (because j<i) without the diagonal, multiply by 2,
    #         % sum.
    #         % transpose it to have it as a column vector
             γ = sum(2 .* S, dims=1);

             return vec(γ)

end