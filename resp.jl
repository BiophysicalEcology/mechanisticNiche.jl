


#### RESPIRATION

include("aux_functions.jl")

function resp(;
    Tskin = 23, 
    Ww_g = 0.04, 
    Tc = 25, 
    Qmet = 0.01241022, 
    EEf = 20, # EXTREF extraction efficiency
    pant = 1, 
    RQ = 0.8, 
    Ta::Real = 20, 
    rh::Real = 50, 
    bp::Real = 101325, 
    O2gas = 20.95, 
    CO2gas = 0.03, 
    N2gas = 79.02, 
    ΔT = 0.1 # DELTAR: temperature difference (ºC) between expired and inspired air # currently not in use?
    )

    Tair = Ta

    rpctO2 = 0.2095
    rpctCO2 = 3e-04
    rpctN2 = 0.7902
    
    
    if (rpctO2 != O2gas/100)
        pctO2 = O2gas/100
    else 
        pctO2 = rpctO2
    end

    if (rpctCO2 != CO2gas/100)
        pctCO2 = CO2gas/100
    else
        pctCO2 = rpctCO2
    end

    if (rpctN2 != N2gas/100)
        pctN2 = N2gas/100
    else
        pctN2 = rpctN2
    end

    TOTgas = pctO2 + pctN2 + pctCO2
    if (TOTgas > 1)
        pctO2 = 1 - (pctN2 + pctCO2)
    elseif (TOTgas < 1) 
        pctO2 = 1 - (pctN2 + pctCO2)
    end

    RGC = 8314.46 # R Gas Constant (RGC): Universal Gas constant (PA - L)/(MOL - K)
    wb = 0
    dp = 999

    pO2 = bp * pctO2
    refpO2 = 101325 * rpctO2
    
    XCALC = Tskin

    if (Tc > 50)
        XCALC = 50
    elseif (Tc < 0)
        XCALC = 0.01
    end

    
    O2stp = (1 / 3.6e6) * (Qmet / 0.0056)
    
    Tlung = XCALC

    VO2 = (O2stp * pO2 / 273.15) * ((Tlung + 273.15) / pO2) # VOLUME OF O2 consumed AT TLUNG
    O2mol = bp * VO2 / (RGC * (Tlung + 273.15))
    O2mol_ext = O2mol / (EEf / 100)
    N2mol_ext = O2mol_ext * (pctN2 / pctO2)
    
    Vair = VO2 / pctO2
    
    VCO2 = pctCO2 * Vair
    CO2mol = bp * VCO2 / (RGC * (Tlung + 273.15))
    
    AIRmol = (O2mol_ext + N2mol_ext + CO2mol) * pant
    AIRvol = (AIRmol * RGC * 273.15 / 101325)
    
    wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol = AIRmol * (esat * (rh / 100)) / (bp - esat * (rh / 100))
    O2mol1 = O2mol_ext - O2mol # remove consumed oxygen from the total
    N2mol1 = N2mol_ext
    CO2mol1 = RQ * O2mol + CO2mol
    AIRmol1 = (O2mol1 + N2mol1 + CO2mol1) * pant
    
    db = XCALC
    rh_exit = 100
    
    wair = wetair(db=db, wb=wb, rh=rh_exit, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol1 = AIRmol1 * (esat / (bp - esat))
    EVPmol = Wmol1 - Wmol
    
    Gevap = EVPmol * 18 # grams lost evaporated through breathing
    
    Kgevap = Gevap / 1000
    Hvpr = 2501200 - 2378.7 * Tlung
    Qresp = Hvpr * Kgevap

    return (Qresp = Qresp, GEVAP = Gevap, AIRML1 = AIRmol, 
        AIRML2 = AIRmol1, WMOL1 = Wmol, WMOL2 = Wmol1, O2MOL1 = O2mol_ext, 
        O2MOL2 = O2mol1, CO2MOL1 = CO2mol, CO2MOL2 = CO2mol1)

end



function resp(o::fixedParams,
    e::fixedEnvParams,
    v::envVars,
    s::stateVariables,
    Qmet = 0.01241022)

    @unpack Ww_g, EEf, RQ, pant = o
    @unpack bp, O2gas, CO2gas, N2gas = e
    @unpack Ta, rh = v
    @unpack Tc, Tskin = s
    

    Tair = Ta[1]
    rh = rh[1]

    rpctO2 = 0.2095
    rpctCO2 = 3e-04
    rpctN2 = 0.7902
    
    
    if (rpctO2 != O2gas/100)
        pctO2 = O2gas/100
    else 
        pctO2 = rpctO2
    end

    if (rpctCO2 != CO2gas/100)
        pctCO2 = CO2gas/100
    else
        pctCO2 = rpctCO2
    end

    if (rpctN2 != N2gas/100)
        pctN2 = N2gas/100
    else
        pctN2 = rpctN2
    end

    TOTgas = pctO2 + pctN2 + pctCO2
    if (TOTgas > 1)
        pctO2 = 1 - (pctN2 + pctCO2)
    elseif (TOTgas < 1) 
        pctO2 = 1 - (pctN2 + pctCO2)
    end

    RGC = 8314.46 # R Gas Constant (RGC): Universal Gas constant (PA - L)/(MOL - K)
    wb = 0
    dp = 999

    pO2 = bp * pctO2
    refpO2 = 101325 * rpctO2
    mass_g = Ww_g
    
    XCALC = Tc

    if (Tc > 50)
        XCALC = 50
    elseif (Tc < 0)
        XCALC = 0.01
    end

    
    O2stp = (1 / 3.6e6) * (Qmet / 0.0056)
    
    Tlung = XCALC

    VO2 = (O2stp * pO2 / 273.15) * ((Tlung + 273.15) / pO2) # VOLUME OF O2 consumed AT TLUNG
    O2mol = bp * VO2 / (RGC * (Tlung + 273.15))
    O2mol_ext = O2mol / (EEf / 100)
    N2mol_ext = O2mol_ext * (pctN2 / pctO2)
    
    Vair = VO2 / pctO2
    
    VCO2 = pctCO2 * Vair
    CO2mol = bp * VCO2 / (RGC * (Tlung + 273.15))
    
    AIRmol = (O2mol_ext + N2mol_ext + CO2mol) * pant
    AIRvol = (AIRmol * RGC * 273.15 / 101325)
    
    wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol = AIRmol * (esat * (rh / 100)) / (bp - esat * (rh / 100))
    O2mol1 = O2mol_ext - O2mol # remove consumed oxygen from the total
    N2mol1 = N2mol_ext
    CO2mol1 = RQ * O2mol + CO2mol
    AIRmol1 = (O2mol1 + N2mol1 + CO2mol1) * pant
    
    db = XCALC
    rh_exit = 100
    
    wair = wetair(db=db, wb=wb, rh=rh_exit, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol1 = AIRmol1 * (esat / (bp - esat))
    EVPmol = Wmol1 - Wmol
    
    Gevap = EVPmol * 18 # grams lost evaporated through breathing
    
    Kgevap = Gevap / 1000
    Hvpr = 2501200 - 2378.7 * Tlung
    Qresp = Hvpr * Kgevap

    return (Qresp = Qresp, GEVAP = Gevap, AIRML1 = AIRmol, 
        AIRML2 = AIRmol1, WMOL1 = Wmol, WMOL2 = Wmol1, O2MOL1 = O2mol_ext, 
        O2MOL2 = O2mol1, CO2MOL1 = CO2mol, CO2MOL2 = CO2mol1)

end
    
   

# for zero solver
function resp(Tx, 
    o::fixedParams,
    e::fixedEnvParams,
    v::envVars,
    Qmet::Real)

    Tskin = Tx
    Tc = Tx

    @unpack Ww_g, EEf, RQ, pant = o
    @unpack bp, O2gas, CO2gas, N2gas = e
    @unpack Ta, rh = v
    

    Tair = Ta = Ta[1]
    rh = rh[1]


    rpctO2 = 0.2095
    rpctCO2 = 3e-04
    rpctN2 = 0.7902
    
    
    if (rpctO2 != O2gas/100)
        pctO2 = O2gas/100
    else 
        pctO2 = rpctO2
    end

    if (rpctCO2 != CO2gas/100)
        pctCO2 = CO2gas/100
    else
        pctCO2 = rpctCO2
    end

    if (rpctN2 != N2gas/100)
        pctN2 = N2gas/100
    else
        pctN2 = rpctN2
    end

    TOTgas = pctO2 + pctN2 + pctCO2
    if (TOTgas > 1)
        pctO2 = 1 - (pctN2 + pctCO2)
    elseif (TOTgas < 1) 
        pctO2 = 1 - (pctN2 + pctCO2)
    end

    RGC = 8314.46 # R Gas Constant (RGC): Universal Gas constant (PA - L)/(MOL - K)
    wb = 0
    dp = 999

    pO2 = bp * pctO2
    refpO2 = 101325 * rpctO2
    mass_g = Ww_g
    
    XCALC = Tx

    if (Tx > 50)
        XCALC = 50
    elseif (Tx < 0)
        XCALC = 0.01
    end

    
    O2stp = (1 / 3.6e6) * (Qmet / 0.0056)
    
    Tlung = XCALC

    VO2 = (O2stp * pO2 / 273.15) * ((Tlung + 273.15) / pO2) # VOLUME OF O2 consumed AT TLUNG
    O2mol = bp * VO2 / (RGC * (Tlung + 273.15))
    O2mol_ext = O2mol / (EEf / 100)
    N2mol_ext = O2mol_ext * (pctN2 / pctO2)
    
    Vair = VO2 / pctO2
    
    VCO2 = pctCO2 * Vair
    CO2mol = bp * VCO2 / (RGC * (Tlung + 273.15))
    
    AIRmol = (O2mol_ext + N2mol_ext + CO2mol) * pant
    AIRvol = (AIRmol * RGC * 273.15 / 101325)
    
    wair = wetair(db=Ta, wb=wb, rh=rh, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol = AIRmol * (esat * (rh / 100)) / (bp - esat * (rh / 100))
    O2mol1 = O2mol_ext - O2mol # remove consumed oxygen from the total
    N2mol1 = N2mol_ext
    CO2mol1 = RQ * O2mol + CO2mol
    AIRmol1 = (O2mol1 + N2mol1 + CO2mol1) * pant
    
    db = XCALC
    rh_exit = 100
    
    wair = wetair(db=db, wb=wb, rh=rh_exit, dp=dp, bp=bp)
    esat = wair.esat
    
    Wmol1 = AIRmol1 * (esat / (bp - esat))
    EVPmol = Wmol1 - Wmol
    
    Gevap = EVPmol * 18 # grams lost evaporated through breathing
    
    Kgevap = Gevap / 1000
    Hvpr = 2501200 - 2378.7 * Tlung
    Qresp = Hvpr * Kgevap

    return (Qresp = Qresp, GEVAP = Gevap, AIRML1 = AIRmol, 
        AIRML2 = AIRmol1, WMOL1 = Wmol, WMOL2 = Wmol1, O2MOL1 = O2mol_ext, 
        O2MOL2 = O2mol1, CO2MOL1 = CO2mol, CO2MOL2 = CO2mol1)

end


# Tskin = 23
# mass = 0.04
# Tc = 25
# Qmet = 0.01241022
# EEf = 20
# pant = 1 
# RQ = 0.8 
# Ta = 20 
# rh = 50 
# bp = 101325
# O2gas = 20.95
# CO2gas = 0.03 
# N2gas = 79.02 
# ΔT = 0.1 

# resp(Tskin = Tskin, 
#     mass = mass, 
#     Tc = Tc, 
#     Qmet = Qmet, 
#     EEf = EEf, # EXTREF extraction efficiency
#     pant = pant, 
#     RQ = RQ, 
#     Ta = Ta, 
#     rh = rh, 
#     bp = bp, 
#     O2gas = O2gas, 
#     CO2gas = CO2gas, 
#     N2gas = N2gas, 
#     ΔT = ΔT)





    