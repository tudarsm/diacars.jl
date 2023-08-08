struct MolParams
    WE::Float64      
    WX::Float64      
    WY::Float64      
    WZ::Float64      
    BE::Float64      
    ALPHAE::Float64  
    DE::Float64      
    BETAE::Float64   
    GAME::Float64    
    DELTE::Float64   
    H0::Float64      
    HE::Float64      
    RE::Float64      
    GAM::Float64     
    DGAMDR::Float64  
    GNE::Float64     
    GNO::Float64     
    DFSPHI::Float64  
    AG::Float64      
    CHINR::Float64   
    AC1::Float64     
    CMASS::Float64   
    a1::Float64      
    p2p1_iso::Float64
    p2p1_ani::Float64
end

function getMolecularParameters(species)

    if species == "N2"
        MolParam=MolParams(
            0.235854024E+04,
            0.14305770E+02,
            -0.50668000E-02,
            -0.10950000E-03,
            0.19982600E+01,
            0.17303500E-01,
            0.57740000E-05,
            0.15500000E-07,
            -0.31536099E-04,
            0,
            0.30000000E-11,
            0.18000000E-11,
            0.20743101E+01,
            0.47640000E+01,
            0.74000000E+01,
            0.60000000E+01,
            0.30000000E+01,
            0,
            0.11660000E+01,
            0.85000000E+01/c.NA*c.molvol,   # convert here to chinr per molecule
            0.24158000E-11,
            0.28013000E+02,
            -2.7,
            0.31,
            0.57
        )
    elseif species == "O2"
        MolParam=MolParams(
            0.15800879E+04,
            0.11898533E+02,
            0.30488333E-01,
            0,
            0.14456077E+01,
            0.15879232E-01,
            0.48436250E-05,
            -0.23200000E-08,
            0.33648000E-04,
            0.38000000E-09,
            0.28000000E-11,
            0.,
            0.22890000E+01,
            0.74970000E+01,
            0.74000000E+01,
            0.00000000E+01,
            0.10000000E+01,
            0.50000000E-01,
            0.17350000E+01,
            7.85380000E+00/c.NA*c.molvol,   # convert here to chinr per molecule
            0.25300000E-11,
            0.32000000E+02,
            -3,                # a1, p2p1 from Marrocco, used for calculation of the HW factors DOI 10.1002/jrs.2201
            0.96,
            1.22,
        )
    elseif species == "CO"
        println("CO MODEL NOT VALIDATED")
        MolParam=MolParams(
            0.21698135E+04,
            0.13288310E+02,
            0.10511000E-01,
            0.57400001E-04,
            0.19312809E+01,
            0.17504411E-01,
            0.61206300E-05,
            -0.11500000E-08,
            0.54870000E-06,
            0.18000000E-09,
            0.54765001E-11,
            -0.17300000E-12,
            0.21322169E+01,
            0.363900000E+0,
            0.93000000E+01,
            0.10000000E+01,
            0.10000000E+01,
            0.50000001E-01,
            0.15500000E+01,
            1.23000000E+01/c.NA*c.molvol,
            0.24038100E-11,
            0.28010000E+02,
            -1.,                 # choosing -1 here to force hermann wallis factors to unity. hermann wallis factors are in fact close to 1 for CO, see https://www.sciencedirect.com/science/article/pii/S0022407301000140#BIB17
            0.,
            0.,     
        )
    else
        error("Species $species not defined!")
    end

    return MolParam

end