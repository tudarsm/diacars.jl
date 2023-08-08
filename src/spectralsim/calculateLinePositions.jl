struct Transitions
    O
    Q
    S
end

function calculateLinePositions(T)
    # % for every specimen, calculate the line position based on the CARS
    # % selection rules:
    # % DV = 1 (only ro-vibrational implemented)
    # % DJ = -2,0,2 (O,Q,S branch)
    # % Pure Rotational not yet implemented! This can be easily done if another
    # % case is introduced (let's say, 'ROT'). However, this needs some adaptions
    # % to the later calculations of the amplitudes of Chi(3). In particular, the
    # % Placzek-Teller-coefficients and Hermann-Wallis-Factors have to be
    # % included.

#     % Q branch: DJ = 0
Q = diff(T, dims=2)
#     s.transitions.(specimen{:}).Q = diff(s.mol.(specimen{:}).T,1,2);
#     % O branch: DJ = -2
#     % initiallize with NaNs to have same size arrays with the index still
#     % representing ground state J
O = fill!(similar(Q),NaN)
O[3:end,:] = T[1:end-2,2:end]-T[3:end,1:end-1]
#     s.transitions.(specimen{:}).O = NaN(size(s.transitions.(specimen{:}).Q));
#     % subtract a shifted array to account for DJ=-2.
#     s.transitions.(specimen{:}).O(3:end,:) = s.mol.(specimen{:}).T(1:end-2,2:end)-s.mol.(specimen{:}).T(3:end,1:end-1);
#     % S branch: DJ = +2
#     % same as for O, shift in the other direction
#     s.transitions.(specimen{:}).S = NaN(size(s.transitions.(specimen{:}).Q));
S = fill!(similar(Q),NaN)
#     s.transitions.(specimen{:}).S(1:end-2,:) = s.mol.(specimen{:}).T(3:end,2:end)-s.mol.(specimen{:}).T(1:end-2,1:end-1);
S[1:end-2,:] = T[3:end,2:end]-T[1:end-2,1:end-1]

transitions = Transitions(O,Q,S)

return transitions

end