

# for new custom types


function fun_nmr(Tx,
    o::fixedParams=anim,
    e::fixedEnvParams=env,
    v::envVars=vars,
    s::stateVariables=state,
    QSOLAR=Qsolar,
    QIRIN=Qrin
)

    @unpack Ww_g, shape, RINSUL, M1, M2, M3, EEf, pant, RQ, rho_body = o
    @unpack ψbody, pct_wet, pct_eyes, FAsky, FAsub, FAobj, εa, postur = o
    @unpack εsub, εsky, fluid, O2gas, CO2gas, N2gas, bp, elev = e
    @unpack Ta, Tsky, Tsub, TGRD, vel, rh, ksub = v

    # rename variables
    GEOMETRY = shape
    FLSHCOND = rho_body
    FLTYPE = fluid

    geom = geometry(o)
    VOL = geom.VOL
    R1 = geom.R1

    # set silhouette area
    if postur == 1
        ASIL = geom.ASILN
    elseif postur == 0
        ASIL = (geom.ASILN + geom.ASILP) / 2
    end

    X = Tx

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
        QMETAB = qmet(Tx, o)
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
            RESPout = resp(Tx, o, e, v, QMETAB)
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


    CONV = convection(Tx, o, e, v)
    QCONV = CONV.QCONV
    HD = CONV.HD


    RESPout = resp(Tx, o, e, v, QMETAB)
    QRESP = RESPout.Qresp
    GEVAP = RESPout.GEVAP

    QSEVAP = skinevap(Tx, o, e, v, GEVAP, HD).QSEVAP


    QIROUT = qrout(Tx, o)

    QCOND = qcond(Tx, o, v, 0.025)

    QIN = QSOLAR + QIRIN + QMETAB
    QOUT = QRESP + QSEVAP + QIROUT + QCONV + QCOND
    #C     FINDING THE DEVIATION FROM ZERO IN GUESSING THE SOLUTION
    ENB = QIN - QOUT
    FUN = ENB

    return FUN
end