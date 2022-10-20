
cd("/Users/usuario/Desktop/biophys/0_scripts/julia/simulation")

using CSV, DataFrames
using Roots
using Plots

include("../aux_functions.jl")
include("../qsolar.jl")
include("../qrin.jl")
include("../qmet.jl")
include("../qrout.jl")
include("../resp.jl")
include("../skinevap.jl")
include("../qcond.jl")
include("../convection.jl")
include("../geometry.jl")
include("../temp_heterogen.jl")

function energy_bal(Tskin = Tskin,
    Ta = Ta,
    Tsub = Tsub,
    RH = RH,
    ksub = ksub,
    VEL = VEL,
    Atot = geom.AREA,
    Av = geom.AV,
    Al = geom.AL,
    Aeff = geom.AEFF,
    FAsky = 0.4,
    FAsub = 0.4,
    εa = 0.95,
    qsolar = Qsolar,
    qrin = Qrin,
    # qmet = Qmet,
    AMASS = AMASS,
    M_1 = 0.013,
    M_2 = 0.8,
    M_3 = 0.038
)
    
    # Qin
    Qmet = qmet(Tskin = Tskin, mass = AMASS, M1 = M_1, M2 = M_2, M3 = M_3)

    # Qout
    Qrout = qrout(Tskin = Tskin, Atot = Atot, Av = Av, FAsky = FAsky, FAsub = FAsub, εa = εa)
    Qcond = qcond(Tskin = Tskin, Tsub = Tsub, ksub = ksub, Av=Av)
    Qconv = convection(Tskin = Tskin, TA = Ta, VEL = VEL, ATOT = Atot, AV = Av, AL=Al).QCONV
    Qsevap = skinevap(Tskin = Tskin, Ta = Ta, rh = RH, Aeff=Aeff, Atot = Atot, vel=VEL).QSEVAP
    Qresp = resp(Tskin = Tskin, Ta = Ta, Qmet = Qmet, mass = AMASS, pant=0).Qresp

    (Qsolar + Qrin + Qmet) - (Qrout + Qcond + Qconv + Qsevap + Qresp)
end


# Read environmental data

metout = CSV.read("metout.csv", DataFrame, normalizenames=true)
soil = CSV.read("soil.csv", DataFrame, normalizenames=true)
tcond = CSV.read("tcond.csv", DataFrame, normalizenames=true)


# get required inputs
TAs = metout.TALOC
TGRDs = soil.DEP1
TSKYs = metout.TSKYC
VELs = metout.VLOC
RHs = metout.RHLOC
QSOLRs = metout.SOLR
Zs = metout.ZEN
K_subs = tcond.CD1


tbs = []

for i in 1:length(TAs)

    global Ta = TAs[i]
    global Tsub = TGRDs[i]
    global Tsky = TSKYs[i]
    global VEL = VELs[i]
    global RH = RHs[i]
    global QSOLR = QSOLRs[i]
    global Z = Zs[i]
    global ksub = K_subs[i]


    # geometry
    # geom = geometry(AMASS = 0.03)# check why not working
    global AMASS = 0.03
    global PTCOND = 0.3
    global geom = geometry(AMASS = AMASS, PTCOND = PTCOND)
    
    # QSOLAR

    ABSAN = 0.85
    ABSSB = 1 - 0.15
    FATOSK = FAsky = 0.4
    FATOSB = FAsub = 0.4
    FATOBJ = FAobj = 0
    ZEN = Z
    PDIF = 0.1
    SHADE = 0
    global POSTUR = postur = 0

    global Qsolar = qsolar(Atot=geom.AREA,
        Asil=geom.ASILN,
        Av=geom.AV,
        At=geom.AT,
        αa=ABSAN, # ABSAN
        αsub=ABSSB, # ABSSB
        FAsky=FATOSK,
        FAsub=FATOSB,
        FAobj=FATOBJ,
        zen=ZEN,
        QSOLR=QSOLR, # incoming solar radiation?
        Pdif=PDIF, # proportion of solar radiation that is diffuse (fractional, 0-1)
        shade=SHADE,
        postur=POSTUR).Qsolar


    # QRIN

    εa = 0.95
    εsub = 0.95
    εsky = 0.8
    
    global Qrin = qrin(Atot=geom.AREA,
        Av=geom.AV,
        At=geom.AT,
        FAsky=FAsky,
        FAsub=FAsub,
        FAobj=FAobj,
        εa=εa,
        εsub=εsub,
        εsky=εsky,
        Tsky=Tsky,
        Tsub=Tsub)
    

    # QMET

    # global M_1 = 0.013
    #M_2 = 0.8
    #M_3 = 0.038
    # M_1 = 0.0 # turn off metabolic heat
    # M_2 = 0.8
    # M_3 = 0.038

    # global Qmet = qmet(Tskin = Ta, mass = AMASS, M1 = M_1, M2 = M_2, M3 = M_3)

    tb_new = find_zero(energy_bal, (-50, 100), Bisection())
    # tb_new = fzero(energy_bal, (-50, 100), Roots.Brent())
    # tb_new = fzero(energy_bal, -50, 100)
    
    push!(tbs, tb_new)

end

maximum(tbs)
minimum(tbs)

plot(tbs)
vline!([152])
hline!([20])



