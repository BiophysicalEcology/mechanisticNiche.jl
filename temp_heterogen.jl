

# from core temperature (Tc) calculate skin (Tskin) and lung (Tlung) temperature

# only running for lizard shape

function temp_heterogen(;
    TC,
    GEOMETRY = 3,
    VOL,
    R1,
    RINSUL,
    FLSHCOND,
    QMETAB,
    QRESP
)

    QGENET = QMETAB - QRESP
    GN = QGENET / VOL
    if GEOMETRY == 0
        TSKIN = TC - GN * R^2 / (2 * FLSHCOND)
        RFLESH = R
        TLUNG = TC
    elseif GEOMETRY == 1
        RFLESH = R1 - RINSUL
        TSKIN = TC - GN * RFLESH^2 / (4 * FLSHCOND)
        TLUNG = (GN * RFLESH^2) / (8 * FLSHCOND) + TSKIN
    elseif GEOMETRY == 2
        A = ASEMAJR
        B = BSEMINR
        C = CSEMINR
        ASQ = A^2
        BSQ = B^2
        CSQ = C^2
        TSKIN = TC - (GN / (2 * FLSHCOND)) * ((ASQ * BSQ * CSQ) / (ASQ * BSQ + ASQ * CSQ + BSQ * CSQ))
        TLUNG = (GN / (4 * FLSHCOND)) * ((ASQ * BSQ * CSQ) / (ASQ * BSQ + ASQ * CSQ + BSQ * CSQ)) + TSKIN
    elseif GEOMETRY == 4
        RFLESH = R1 - RINSUL
        RSKIN = R1
        S1 = (QGENET / (4 * PI * FLSHCOND)) * ((RFLESH - RSKIN) / (RFLESH * RSKIN))
        TSKIN = TC - (GN * RFLESH^2) / (6 * FLSHCOND) + S1
        TLUNG = (GN * RFLESH^2) / (12 * FLSHCOND) + TSKIN
    elseif GEOMETRY == 3 || GEOMETRY == 5
        RFLESH = R1 - RINSUL
        TSKIN = TC - GN * RFLESH^2 / (4 * FLSHCOND)
        TLUNG = (GN * RFLESH^2) / (8 * FLSHCOND) + TSKIN
    end

    if TLUNG > TC
        TLUNG = TC
    end

    return (Tskin = TSKIN, Tlung = TLUNG)
end


# include("geometry.jl")
# include("qmet.jl")
# include("resp.jl")

# geom = geometry()

# Tc = TC = 25
# AMASS = 40
# M_1 = 0.013
# M_2 = 0.8
# M_3 = 0.038
# Qmet = qmet(Tskin = Tc, mass = AMASS, M1 = M_1, M2 = M_2, M3 = M_3)
# Qresp = resp(Tskin = Tc).Qresp

# temp_heterogen(TC = Tc,
#     GEOMETRY = 3, 
#     R1 = geom.R1,
#     VOL = geom.VOL,
#     RINSUL = 0,
#     FLSHCOND = 0.5,
#     QMETAB = Qmet,
#     QRESP = Qresp)