struct Amplitudes
    O
    Q
    S
end

function chiamplitudes(MolParam,Δ)
# % this function computes the amplitudes of the resonant susceptibility
# % note: at this point, this is without the number density N, which will be
# % treated later
# for specimen = s.species
#     % constants
#     zeta   = s.mol.(specimen{:}).AC1;
#     vp1 = s.V(2:end);   % this is v+1, i.e. ground state +1
#     xi = s.mol.(specimen{:}).AG;

zeta = MolParam.AC1
vp1 = V[2:end]
xi = MolParam.AG

#     % Placzek-Teller coefficients
#     bjjO = s.J.*(s.J-1)./((2*s.J-1).*(2*s.J+1));
#     bjjQ = s.J.*(s.J+1)./((2*s.J-1).*(2*s.J+3));
#     bjjS = (s.J+1).*(s.J+2)./((2*s.J+3).*(2*s.J+1));

bjjO = @. J*(J-1)/((2*J-1)*(2*J+1))
bjjQ = @. J*(J+1)/((2*J-1)*(2*J+3))
bjjS = @. (J+1)*(J+2)/((2*J+3)*(2*J+1))

# println(bjjO)
# println(bjjQ)
# println(bjjS)

#     % Herrmann-Wallis factors
eta = 2*MolParam.BE/MolParam.WE;    # eta in Marrocco, gamma in Klemperer
#     % for O and S-branch
eta_os = 4*MolParam.BE*MolParam.GAM/(MolParam.WE*MolParam.DGAMDR*MolParam.RE); # this is from the left central part in Buckingham, p47. This corresponds to 4*BE/WE*(alpha||-alpha_|_)_e/(alpha||-alpha_|_)'_e
FS = @. (1-eta_os.*(2*J+3)).^2;    #% buckingham, p47
FO = @. (1+eta_os*(2*J-1)).^2;    #% buckingham, p47

# FOR NOW: HERMANN-WALLIS ONLY TB model
a1 = MolParam.a1;  #% see Marrocco paper, Table1
p2p1_iso = MolParam.p2p1_iso;    #% see Marrocco paper, Table1
p2p1_ani = MolParam.p2p1_ani;    #% see Marrocco paper, Table1
FQ_iso = @. 1-(3*(a1+1)/2-4*p2p1_iso)*eta^2*J.*(J+1);
FQ_ani = @. 1-(3*(a1+1)/2-4*p2p1_ani)*eta^2*J.*(J+1);

#     % Hermann-Wallis-Factors depending on chosen model
#     switch s.HWFactors
#         case 'JK'
#             FQ = 1-3*eta^2*s.J.*(s.J+1)/2; % klemperer, below equation (14)
#         case 'LBY'
#             a1 = s.mol.(specimen{:}).a1; % short handle to a1 constant
#             FQ = (1-3*eta^2*(a1+1)*s.J.*(s.J+1)/4).^2;
#         case 'TB'
#             a1 = s.mol.(specimen{:}).a1;  % see Marrocco paper, Table1
#             p2p1_iso = s.mol.(specimen{:}).p2p1_iso;    % see Marrocco paper, Table1
#             p2p1_ani = s.mol.(specimen{:}).p2p1_ani;    % see Marrocco paper, Table1
#             FQ_iso = 1-(3*(a1+1)/2-4*p2p1_iso)*eta^2*s.J.*(s.J+1);
#             FQ_ani = 1-(3*(a1+1)/2-4*p2p1_ani)*eta^2*s.J.*(s.J+1);
#         otherwise
#             error('Hermann-Wallis-Factors should either be JK, LBY or TB')
#     end

#     % O-branch
Opar = @. 2/15*bjjO.*FO*xi^2;
Operp = @. 3/4*Opar;

#     % Q-branch
#     switch s.HWFactors
#         case 'TB'
Qpar = @. (1+4/45*xi^2*bjjQ).*FQ_iso;
Qperp = @. 1/15*xi^2*bjjQ.*FQ_ani;
#         otherwise
#             Qpar = (1+4/45*xi^2*bjjQ).*FQ;
#             Qperp = 1/15*xi^2*bjjQ.*FQ;
#     end
#     % S-branch
Spar = @. 2/15*bjjS.*FS*xi^2;
Sperp = @. 3/4*Spar;

#     % compute the amplitudes
#     s.chiamp.(specimen{:}).O = 1e18 * s.delta.(specimen{:}).O .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Opar'+sind(s.theta)*sind(s.phi)*Operp');
#     s.chiamp.(specimen{:}).Q = 1e18 * s.delta.(specimen{:}).Q .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Qpar'+sind(s.theta)*sind(s.phi)*Qperp');
#     s.chiamp.(specimen{:}).S = 1e18 * s.delta.(specimen{:}).S .* vp1 * zeta^2/(4*pi*sconst('c')).*(cosd(s.theta)*cosd(s.phi)*Spar'+sind(s.theta)*sind(s.phi)*Sperp');
# end

# so far only parallel polarization
# println(size(Opar))
# println(size(Δ.Q))
aO = @. 1e18 * Δ.O * vp1' * zeta^2 / (4*pi*c.c) * Opar
aQ = @. 1e18 * Δ.Q * vp1' * zeta^2 / (4*pi*c.c) * Qpar
aS = @. 1e18 * Δ.S * vp1' * zeta^2 / (4*pi*c.c) * Spar

return Amplitudes(aO,aQ,aS)

end