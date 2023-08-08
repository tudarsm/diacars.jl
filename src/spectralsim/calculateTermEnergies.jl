
function calculateTermEnergies(MolParam)
#     % Anharmonic oscillator, see following book, which should be available
#     % on google books. The equations are right below figure 47.
#     %     @book{herzberg2013molecular,
#     %   title={Molecular spectra and molecular structure},
#     %   author={Herzberg, Gerhard},
#     %   volume={1},
#     %   year={2013},
#     %   publisher={Read Books Ltd}
#     % }

#     % kept upper case W for omega for consistency, this is equation
#     % (III,77) in herzberg
W0 = MolParam.WE-MolParam.WX+3/4*MolParam.WY
W0X0 = MolParam.WX-3/2*MolParam.WY;
W0Y0 = MolParam.WY;

#     % set the nuclear spin as degeneracy for even and odd Js
#     % ATTENTION: J of course starts counting at zero, matlab at one
#     % so 1,3,5 ... corresponds to the EVEN Js (0,2,4,...)
GJ=zeros(Float64,1,length(J))
GJ[1:2:end].=MolParam.GNE
GJ[2:2:end].=MolParam.GNO

#     % this is equation (III,76)
#     MolParam.G0V = MolParam.W0.*(s.V) - MolParam.W0X0*(s.V).^2 + MolParam.W0Y0*(s.V).^3;
G0V = W0.*V - W0X0.*(V.^2) + W0Y0.*V.^3
#     % Vibrating rotator
#     % See Herzberg, see equations (III,124-126)
#     % to make it a little easier to read, some shortcuts to (v+1/2) and
#     % J(J+1). Also, this is vectorized, make sure that v and J are in
#     % different dimensions.
#     % again, keep the uppercase lettering for easier readability.
     VPH = transpose(V.+0.5);
     JJ = (J.*(J.+1))
     Bv = MolParam.BE.+VPH.*(-MolParam.ALPHAE.+VPH.*MolParam.GAME)
     Dv = MolParam.DE.+VPH.*(MolParam.BETAE.+VPH.*MolParam.DELTE)
     Hv = MolParam.H0.+VPH.*MolParam.HE
     FJ = Bv.*JJ.-Dv.*JJ.^2 .+ Hv.*JJ.^3
#     % Term values for the rotating vibrator (III-121), higher orders of vph
#     % available in CARS.MOL, use all. Note: WE corresponds to omega_e*x_e
#     % and so forth.
     T = FJ .+ MolParam.WE.*VPH .- MolParam.WX.*VPH.^2 .+ MolParam.WY.*VPH.^3 .+ MolParam.WZ.*VPH.^4

     return T,FJ,GJ,G0V
end