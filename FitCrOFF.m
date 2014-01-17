function MRS_struct = FitCrOFF(MRS_struct)
    %This fits Cr in the OFF spectrum. This is used for initialisation of
    %further fitting for freq/phase correction

    %Determine which is OFF
    switch MRS_struct.p.ONOFForder
        case 'ONfirst'
       ToBeFit=MRS_struct.
        case 'OFFfirst'
            
    end
end