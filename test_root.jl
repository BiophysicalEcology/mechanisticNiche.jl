
include("aux_functions.jl")
include("qsolar.jl")
include("qrin.jl")
include("qmet.jl")
include("qrout.jl")
include("resp.jl")
include("skinevap.jl")
include("qcond.jl")
include("convection.jl")
include("geometry.jl")
include("temp_heterogen.jl")


# geometry
geom = geometry()

# QSOLAR

ABSAN = 0.85
ABSSB = 0.2 
FATOSK = FAsky = 0.4
FATOSB = FAsub = 0.4
FATOBJ = FAobj = 0
ZEN = 20 * pi/180
QSOLR = 1000
PDIF = 0.1
SHADE = 0
POSTUR = 1

Qsolar = qsolar(Atot = geom.AREA, 
    Asil = geom.ASILN, 
    Av = geom.AV, 
    At = geom.AT, 
    αa = ABSAN, # ABSAN
    αsub = ABSSB, # ABSSB
    FAsky = FATOSK, 
    FAsub = FATOSB, 
    FAobj = FATOBJ, 
    zen = ZEN, 
    QSOLR = QSOLR, # incoming solar radiation?
    Pdif = PDIF, # proportion of solar radiation that is diffuse (fractional, 0-1)
    shade = SHADE, 
    postur = POSTUR).Qsolar


# QRIN

εa = 0.95
εsub = 0.95
εsky = 0.8
Tsky = -5
Tsub = 30

Qrin = qrin(Atot = geom.AREA, 
    Av = geom.AV, 
    At = geom.AT, 
    FAsky = FAsky,
    FAsub = FAsub, 
    FAobj = FAobj, 
    εa = εa,
    εsub = εsub,
    εsky = εsky,
    Tsky = Tsky, 
    Tsub = Tsub)


# # QMET

# AMASS = 0.04
# XTRY = 25
# M_1 = 0.013
# M_2 = 0.8
# M_3 = 0.038


# Qmet = qmet(mass = AMASS, 
#     XTRY = XTRY,
#     M1 = M_1, 
#     M2 = M_2, 
#     M3 = M_3)




function energy_bal(Tskin = Tskin,
    Ta = Ta,
    Tsub = Tsub,
    Atot = 0.01325006,
    Av = 0.001325006,
    FAsky = 0.4,
    FAsub = 0.4,
    εa = 0.95,
    qsolar = Qsolar,
    qrin = Qrin,
    AMASS = 0.04,
    M_1 = 0.013,
    M_2 = 0.8,
    M_3 = 0.038
)
    
    # Qin
    Qmet = qmet(Tskin = Tskin, mass = AMASS, M1 = M_1, M2 = M_2, M3 = M_3)

    # Qout
    Qrout = qrout(Tskin = Tskin, Atot = Atot, Av = Av, FAsky = FAsky, FAsub = FAsub, εa = εa)
    Qcond = qcond(Tskin = Tskin, Tsub = Tsub, ksub = 0.5)
    Qconv = convection(Tskin = Tskin, TA = Ta).QCONV
    Qsevap = skinevap(Tskin = Tskin, Ta = Ta, rh = 5).QSEVAP
    Qresp = resp(Tskin = Tskin, Ta = Ta).Qresp

    (Qsolar + Qrin + Qmet) - (Qrout + Qcond + Qconv + Qsevap + Qresp)
end


using Roots

Ta = 20
Tsub = 30

find_zero(energy_bal, (-50, 100), Bisection())


# Checking with Mike's R code
Qsolar = 7.563969
Qrin = 3.645377
Qmet = 0.03966429

Tc = round(find_zero(energy_bal, (-50, 100), Bisection()), digits=4)


temp_heterogen(TC = Tc,
    GEOMETRY = 3, 
    R1 = geom.R1,
    VOL = geom.VOL,
    RINSUL = 0,
    FLSHCOND = 0.5,
    QMETAB = Qmet,
    QRESP = Qresp)