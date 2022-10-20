
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


Ww_g = 40
alpha = 0.85
epsilon = 0.95
rho_body = 1000
fatosk = 0.4
fatosb = 0.4
shape = 3
shape_a = 1
shape_b = 3
shape_c = 2 / 3
custom_shape = [10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743]
pct_cond = 10
pct_touch = 0
postur = 1
k_flesh = 0.5
M_1 = 0.013
M_2 = 0.8
M_3 = 0.038
pct_wet = 0.1
pct_eyes = 0
pct_mouth = 0
psi_body = -700
pantmax = 1
F_O2 = 20
RQ = 0.8
delta_air = 0.1
elev = 0
alpha_sub = 0.2
epsilon_sub = 1
epsilon_sky = 1
pres = 101325
fluid = 0
O2gas = 20.95
CO2gas = 0.03
N2gas = 79.02
K_sub = 0.5
PDIF = 0.1
SHADE = 0
QSOLR = 1000
Z = 20
TA = 20
TGRD = 30
TSUBST = 30
TSKY = -5
VEL = 1
RH = 5


# parameter name translations from input arguments, and value conversions

ZEN = Z / 180 * pi
AMASS = Ww_g / 1000 # animal wet weight (kg)
SKINW = pct_wet / 100 # fractional skin wetness
PEYES = pct_eyes / 100 # fractional of surface area that is wet eyes
PMOUTH = pct_mouth / 100 # fraction of surface area that is wet mouth
SKINT = pct_touch / 100 # fraction of surface area that is touching another individual at the same temperature
PTCOND = pct_cond / 100 # fraction of surface area conducting to the ground
ALT = elev # 'altitude' (m) (technically correct term is elevation)
EMISSK = epsilon_sky # emissivity of the sky (0-1)
EMISSB = epsilon_sub # emissivity of the substrate (0-1)
ABSSB = alpha_sub # solar absorbtivity of the substrate (0-1)
BP = pres # barometric pressure
RELHUM = RH # relative humidity
PSI_BODY = psi_body #
ABSAN = alpha # animal solar absorbtivity
EXTREF = F_O2 # oxygen extraction efficiency (0-1)
GEOMETRY = shape # animal shape
FATOSK = fatosk # radiation configuration factor to sky
FATOSB = fatosb # radiation configuration factor to ground
O2GAS = O2gas # O2 gas concentration (%)
CO2GAS = CO2gas # CO2 gas concentration (%)
N2GAS = N2gas # N2 gas concentration (%)
FLSHCOND = k_flesh # flesh thermal conductivity (W/mK)
ANDENS = rho_body # body density (kg/m3)
PANT = pantmax # panting modifier, >=1 (or 0 if cutting out respiration)
FLTYPE = fluid # air or water? (0 or 1)
EMISAN = epsilon # emissivity of animal skin (0-1)
EMISSB = epsilon_sub # emissivity of substrate (0-1)
EMISSK = epsilon_sky # emissivity of sky (0-1)
SUBTK = K_sub # substrate thermal conductivity (W/mK)
DELTAR = delta_air # temperature difference between inspired and expired air
CUSTOMGEOM = custom_shape # parameters for customised geometry

# unused parameters
FATOBJ = 0 # configuration factor to nearby object of different temp to sky and ground (e.g. warm rock, fire), not yet used
RINSUL = 0 # radius of insulation, not used yet

# shape setup for custom animals
if shape == 3
    shape_a = 1
    shape_b = 1
    shape_c = 4
end
# if(shape == 4){ # frog proportions
#   shape_a = 1
#   shape_b = 1
#   shape_c = 0.5
# }
SHP = [shape_a, shape_b, shape_c]
if SHP[1] == SHP[2] && SHP[2] == SHP[3]
    SHP[3] = SHP[3] - 0.0000001
end

# call GEOM_ecto to get lengths, areas and volume
geom = geometry(AMASS=AMASS,
    GEOMETRY=GEOMETRY,
    SHP=SHP,
    CUSTOMGEOM=CUSTOMGEOM,
    ANDENS=ANDENS,
    SKINW=SKINW,
    SKINT=SKINT,
    RINSUL=RINSUL,
    PTCOND=PTCOND,
    PMOUTH=PMOUTH,
    PANT=PANT)


VOL = geom.VOL
AREA = geom.AREA
ATOT = AREA
AV = geom.AV
AT = geom.AT
AL = geom.AL
ASILN = geom.ASILN
ASILP = geom.ASILP
AEFF = geom.AEFF
R1 = geom.R1
R = geom.R
ASEMAJR = geom.ASEMAJR
BSEMINR = geom.BSEMINR
CSEMINR = geom.CSEMINR

# set silhouette area
if postur == 1
    ASIL = ASILN
end
#   if(postur == 2){
#     ASIL = ASILP
#   }
#   if(postur == 0){
#     ASIL = (ASILN + ASILP) / 2
#   }

# compute solar load
SOLAR = qsolar(Atot=AREA,
    Asil=ASIL,
    Av=AV,
    At=AT,
    αa=ABSAN,
    αsub=ABSSB,
    FAsky=FATOSK,
    FAsub=FATOSB,
    FAobj=FATOBJ,
    zen=ZEN,
    QSOLR=QSOLR,
    Pdif=PDIF,
    shade=SHADE,
    postur=postur)
QSOLAR = SOLAR.Qsolar


# compute infrared radiation in
QIRIN = qrin(Atot=ATOT,
    Av=AV,
    At=AT,
    FAsky=FATOSK,
    FAsub=FATOSB,
    FAobj=FATOBJ,
    εa=EMISAN,
    εsub=Float64(EMISSB),
    εsky=Float64(EMISSK),
    Tsky=TSKY,
    Tsub=TGRD)



# going to fun_ecto

X = 25

#     CONTROL OF BODY TEMPERATURE GUESSES FOR STABILITY PURPOSES
if X > 100
    X = 100
end

TC = X
XTRY = X

#C     GET THE METABOLIC RATE
#C     CHECKING FOR INANIMATE OBJECT
#C      ALIVE, BUT IS IT TOO COLD?
if TC >= 0
    QMETAB = qmet(mass=AMASS,
        Tskin=XTRY,
        M1=M_1,
        M2=M_2,
        M3=M_3)
    # QMETAB = MET.out
else
    #C       TOO COLD, SUPER LOW METABOLISM
    QMETAB = 0.0001
    TC = X
end


#C     GET THE RESPIRATORY WATER LOSS
#C     CHECKING FOR FLUID TYPE
if FLTYPE == 0
    #C      AIR
    #C      CALL FOR RESPIRATORY WATER & ENERGY LOSS
    if QMETAB >= 0
        RESPout = resp(Tskin=XTRY,
            mass=AMASS,
            Tc=TC,
            Qmet=QMETAB,
            EEf=EXTREF,
            pant=PANT,
            RQ=RQ,
            Ta=TA,
            rh=RELHUM,
            bp=BP,
            O2gas=O2GAS,
            CO2gas=CO2GAS,
            N2gas=N2GAS)
        QRESP = RESPout.Qresp
    else
        #C       NEGATIVE METABOLIC RATE. NO PHYSIOLOGICAL MEANING - DEAD.
        QRESP = 0
        QMETAB = 0
    end
end



#C     NET INTERNAL HEAT GENERATION
QGENET = QMETAB - QRESP
#C     NET INTERNAL HEAT GENERATION/UNIT VOLUME. USE FOR ESTIMATING SKIN TEMP.
GN = QGENET / VOL

#C     COMPUTING SURFACE TEMPERATURE AS DICTATED BY GEOMETRY
#C     FIRST SET AVERAGE BODY TEMPERATURE FOR ESTIMATION OF AVEARAGE LUNG TEMPERATURE
if GEOMETRY == 3 || GEOMETRY == 5
    # C      MODEL LIZARD/CUSTOM SHAPE AS CYLINDER
    # C      CYLINDER: FROM P. 270 BIRD, STEWART & LIGHTFOOT. 1960. TRANSPORT PHENOMENA.
    # C      TAVE = (GR ^ 2/(8K)) + TSKIN, WHERE TSKIN = TCORE - GR ^ 2/(4K)
    # C      NOTE:  THESE SHOULD ALL BE SOLVED SIMULTANEOUSLY.  THIS IS AN APPROXIMATION
    # C      USING CYLINDER GEOMETRY. SUBCUTANEOUS FAT IS ALLOWED IN CYLINDER & SPHERE
    # C      CALCULATIONS.
    RFLESH = R1 - RINSUL
    TSKIN = TC - GN * RFLESH^2 / (4 * FLSHCOND)
    #C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
    TLUNG = (GN * RFLESH^2) / (8 * FLSHCOND) + TSKIN
end

#C     LIMITING LUNG TEMPERATURE EXTREMES
if TLUNG > TC
    TLUNG = TC
end


CONV = convection(GEOMETRY=GEOMETRY,
    ATOT=ATOT,
    AV=AV,
    AL=AL,
    AT=AT,
    BP=BP,
    ALT=ALT,
    TA=TA,
    VEL=VEL,
    # FLTYPE = FLTYPE,
    Tskin=TSKIN)
QCONV = CONV.QCONV
HD = CONV.HD


RESPout = resp(Tskin=XTRY,
    mass=AMASS,
    Tc=TC,
    Qmet=QMETAB,
    EEf=EXTREF,
    pant=PANT,
    RQ=RQ,
    Ta=TA,
    rh=RELHUM,
    bp=BP,
    O2gas=O2GAS,
    CO2gas=CO2GAS,
    N2gas=N2GAS)

QRESP = RESPout.Qresp
GEVAP = RESPout.GEVAP

SEVAPout = skinevap(#TC = TC,
    Tskin=TSKIN,
    Gevap=GEVAP,
    ψbody=PSI_BODY,
    skinwet=SKINW,
    Aeff=AEFF,
    Atot=ATOT,
    hd=HD,
    peyes=PEYES,
    Ta=TA,
    rh=RELHUM,
    vel=VEL,
    bp=BP)
QSEVAP = SEVAPout.QSEVAP


QIROUT = qrout(Tskin=TSKIN,
    Atot=ATOT,
    Av=AV,
    # AT = AT,
    FAsky=FATOSK,
    FAsub=FATOSB,
    εa=EMISAN)

QCOND = qcond(Av=AV,
    Tskin=TSKIN,
    Tsub=TSUBST,
    ksub=SUBTK)


QIN = QSOLAR + QIRIN + QMETAB
QOUT = QRESP + QSEVAP + QIROUT + QCONV + QCOND
#C     FINDING THE DEVIATION FROM ZERO IN GUESSING THE SOLUTION
ENB = QIN - QOUT
FUN = ENB


