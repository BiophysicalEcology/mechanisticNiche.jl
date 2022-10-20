

function fun_nmr(X,
    AMASS = AMASS,
    GEOMETRY = GEOMETRY,
    ATOT = AREA,
    AV = AV,
    AT = AT,
    AL = AL,
    VOL = VOL,
    R = R,
    R1 = R1,
    RINSUL = RINSUL,
    ASEMAJR = ASEMAJR,
    BSEMINR = BSEMINR,
    CSEMINR = CSEMINR,
    M_1 = M_1,
    M_2 = M_2,
    M_3 = M_3,
    EXTREF = EXTREF,
    PANT = PANT,
    RQ = RQ,
    FLSHCOND = FLSHCOND,
    PSI_BODY = PSI_BODY,
    SKINW = SKINW,
    AEFF = AEFF,
    PEYES = PEYES,
    FATOSK = FATOSK,
    FATOSB = FATOSB,
    FATOBJ = FATOBJ,
    EMISAN = EMISAN,
    EMISSB = EMISSB,
    EMISSK = EMISSK,
    FLTYPE = FLTYPE,
    TA = TA,
    TSKY = TSKY,
    TSUBST = TSUBST,
    TGRD = TGRD,
    VEL = VEL,
    QSOLAR = QSOLAR,
    QIRIN = QIRIN,
    RELHUM = RELHUM,
    BP = BP,
    ALT = ALT,
    SUBTK = SUBTK,
    O2GAS = O2GAS,
    CO2GAS = CO2GAS,
    N2GAS = N2GAS)

    #     CONTROL OF BODY TEMPERATURE GUESSES FOR STABILITY PURPOSES
    if X > 100
        X = 100
    end

    TC = X
    XTRY = X

    Ww_g = AMASS * 1000

    #C     GET THE METABOLIC RATE
    #C     CHECKING FOR INANIMATE OBJECT
    #C      ALIVE, BUT IS IT TOO COLD?
    if TC >= 0
        QMETAB = qmet(Ww_g=Ww_g,
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
                Ww_g=Ww_g,
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
        Ww_g=Ww_g,
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

    QSEVAP = skinevap(#TC = TC,
        Tskin=TSKIN,
        Gevap=GEVAP,
        ψbody=PSI_BODY,
        pct_wet=SKINW,
        Aeff=AEFF,
        Atot=ATOT,
        hd=HD,
        pct_eyes=PEYES,
        Ta=TA,
        rh=RELHUM,
        vel=VEL,
        bp=BP).QSEVAP
    
        
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

    return FUN
end

