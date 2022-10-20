


## CONVECTION

function convection(;
    Tskin=25,
    GEOMETRY=3,
    ATOT=0.01325006,
    AV=0.001325006,
    AL=0.03419952,
    AT=0,
    BP=101325,
    ALT=0,
    TA=20,
    VEL=1
    # FLTYPE = 0
)



    #C     SETTING THE CHARACTERISTIC DIMENSION FOR NUSSELT-REYNOLDS CORRELATIONS
    D = AL

    BETA = 1 / (TA + 273)
    CP = 1.0057E+3
    G = 9.80665

    #C     CONVECTIVE AREA CALCULATION = TOTAL AREA - VENTRAL AREA IN CONTACT WITH SUBSTRATE
    CONVAR = ATOT - AV - AT


    #C     USING ALTITUDE TO COMPUTE BP (SEE DRYAIR LISTING)
    #C      BP=0.0
    DB = TA

    #C     GET THERMAL PROPERTIES OF DRY AIR AT CURRENT TEMP AND PRESSURE
    dair = dryair(db=DB, bp=BP, alt=ALT)
    PATMOS = dair.patmos
    DENSTY = dair.densty
    VISDYN = dair.visdyn
    VISKIN = dair.viskin
    DIFVPR = dair.difvpr
    THCOND = dair.thcond
    HTOVPR = dair.htovpr
    TCOEFF = dair.tcoeff
    GGROUP = dair.ggroup

    #   #C     CHECKING TO SEE IF THE FLUID IS WATER, NOT AIR
    #   if(FLTYPE == 1){
    #     WATERPROP.out = WATERPROP(TA)
    #     CP = WATERPROP.out$CP
    #     DENSTY = WATERPROP.out$DENSTY
    #     THCOND = WATERPROP.out$THCOND
    #     VISDYN = WATERPROP.out$VISDYN
    #   }

    #C     COMPUTING PRANDTL AND SCHMIDT NUMBERS
    PR = CP * VISDYN / THCOND
    #   if(FLTYPE == 0){
    #     #C      AIR
    SC = VISDYN / (DENSTY * DIFVPR)
    #   }else{
    #     #C      WATER; NO MEANING
    #     SC = 1
    #   }

    #C     SKIN/AIR TEMPERATURE DIFFERENCE
    DELTAT = Tskin - TA
    if DELTAT == 0
        DELTAT = 0.01
    end

    #C     COMPUTING GRASHOF NUMBER
    GR = ((DENSTY^2) * BETA * G * (D^3) * DELTAT) / (VISDYN^2)
    #C     CORRECTING IF NEGATIVE DELTAT
    GR = abs(GR)

    #C     AVOIDING DIVIDE BY ZERO IN FREE VS FORCED RAYLEI
    if VEL <= 0
        VEL = 0.0001
    end

    #C     REYNOLDS NUMBER
    RE = DENSTY * VEL * D / VISDYN

    #C     CHOOSING FREE OR FORCED CONVECTION
    #C     SIGNIFICANT FREE CONVECTION IF GR/RE ^ 2 .GE. 20.0
    #C     KREITH (1965) P. 358
    GRRE2 = GR / (RE^2)

    #C     *********************  FREE CONVECTION  ********************
    #   if(GEOMETRY == 0){
    #     RAYLEI = GR * PR
    #     ANUFRE = 0.55 * RAYLEI ^ 0.25
    #   }

    if GEOMETRY == 1 || GEOMETRY == 3 || GEOMETRY == 5
        #C      FREE CONVECTION OF A CYLINDER
        #C      FROM P.334 KREITH (1965): MC ADAM'S 1954 RECOMMENDED COORDINATES
        RAYLEI = GR * PR
        if RAYLEI < 1.0E-05
            ANUFRE = 0.4
        elseif RAYLEI < 0.1
            ANUFRE = 0.976 * RAYLEI^0.0784
        elseif RAYLEI <= 100
            ANUFRE = 1.1173 * RAYLEI^0.1344
        elseif RAYLEI < 10000.0
            ANUFRE = 0.7455 * RAYLEI^0.2167
        elseif RAYLEI < 1.0e9
            ANUFRE = 0.5168 * RAYLEI^0.2501
        elseif RAYLEI < 1.0e12
            ANUFRE = 0.5168 * RAYLEI^0.2501
        end
    end

    #   if((GEOMETRY == 2) | (GEOMETRY == 4)){
    #     #C      SPHERE FREE CONVECTION
    #     #C      FROM P.413 BIRD ET AL (1960) TRANSPORT PHENOMENA)
    #     RAYLEI = (GR ^ (1 /4)) * (PR ^ (1 / 3))
    #     ANUFRE = 2 + 0.60 * RAYLEI
    #     if(RAYLEI < 200){
    #       #GO TO 20
    #     }else{
    #       message(paste0(RAYLEI, '(GR ^ 0.25) * (PR ^ 0.333) IS TOO LARGE FOR CORREL.'))
    #     }
    #   }

    #C     CALCULATING THE FREE CONVECTION HEAT TRANSFER COEFFICIENT, HC  (NU=HC*D/KAIR)
    HCFREE = (ANUFRE * THCOND) / D

    #C     CALCULATING THE SHERWOOD NUMBER FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SHFREE = ANUFRE * (SC / PR)^(1 / 3)
    #C     CALCULATING THE MASS TRANSFER COEFFICIENT FROM THE SHERWOOD NUMBER
    HDFREE = SHFREE * DIFVPR / D
    #C     CALCULATING THE CONVECTIVE HEAT LOSS AT THE SKIN
    QFREE = HCFREE * CONVAR * (Tskin - TA)

    #C     *******************  FORCED CONVECTION  *********************
    #   if(GEOMETRY == 0){
    #     ANU = 0.102 * RE ^ 0.675 * PR ^ (1 / 3)
    #   }

    #   if(GEOMETRY == 1){
    #     #C      FORCED CONVECTION OF A CYLINDER
    #     #C      ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954)
    #     if(RE < 4){
    #       ANU = 0.891 * RE ^ 0.33
    #     }else{
    #     if(RE < 40){
    #       ANU = 0.821 * RE ^ 0.385
    #     }else{
    #     if(RE < 4000.){
    #       ANU = 0.615 * RE ^ 0.466
    #     }else{
    #     if(RE < 40000.){
    #       ANU = 0.174 * RE ^ 0.618
    #     }else{
    #     if(RE < 400000.){
    #       ANU = 0.0239 * RE ^ 0.805
    #     }}}}}
    #   }

    if GEOMETRY > 1
        #C      FORCED CONVECTION IN SPHERE
        #C       ANU=0.34 * RE ^ 0.24 ! ORIGINAL RELATION
        ANU = 0.35 * RE^0.6 # FROM McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
    end

    #C     FORCED CONVECTION FOR ANIMAL

    # if(GEOMETRY == 4){
    #   #C       ***********************FROG******************************
    #   #C       C.R. TRACY'S LEOPARD FROGS - ECOL. MONOG. 1976 V. 46(3)
    #   #C       CHECKING FOR OUT OF BOUNDS VALUES
    #   if(RE < 80){
    #     message(paste0(' RE, ',RE,',TOO SMALL FOR FROG ANCORR'))
    #   }else{
    #     if (RE > 40000){
    #       message(paste0(' RE, ',RE,',TOO LARGE FOR FROG ANCORR'))
    #     }
    #   }
    #   #C       COMPUTING NUSSELT AND SHERWOOD NUMBERS
    #   if(RE <= 2000){
    #     ANU = 0.88 * RE ^ 0.5
    #   }else{
    #     ANU = 0.258 * RE ^ 0.667
    #   }
    # }

    #C     FORCED CONVECTION FOR ANIMAL
    HCFORC = ANU * THCOND / D # HEAT TRANFER COEFFICIENT
    SHFORC = ANU * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HDFORC = SHFORC * DIFVPR / D # MASS TRANSFER COEFFICIENT
    #C     USING BIRD, STEWART, & LIGHTFOOT'S MIXED CONVECTION FORMULA (P. 445, TRANSPORT PHENOMENA, 2002)
    NUTOTAL = (ANUFRE^3 + ANU^3)^(1 / 3)
    HC = NUTOTAL * (THCOND / D) # MIXED CONVECTION HEAT TRANSFER
    QCONV = HC * CONVAR * (Tskin - TA) # TOTAL CONVECTIVE TRANSFER
    #C     CALCULATING THE SHERWOOD NUMBERs FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SH = NUTOTAL * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HD = SH * DIFVPR / D # MASS TRANSFER COEFFICIENT
    return (CONVAR=CONVAR, QCONV=QCONV, HC=HC, HD=HD, SH=SH, QFREE=QFREE, HCFREE=HCFREE, HCFORC=HCFORC, SHFREE=SHFREE, SHFORC=SHFORC, HDFREE=HDFREE, HDFORC=HDFORC)
end




function convection(o::fixedParams,
    e::fixedEnvParams,
    v::envVars,
    s::stateVariables
)

    @unpack shape = o
    @unpack bp, elev = e
    @unpack Ta, vel = v
    @unpack Tskin = s

    #rename variables: clean the code!
    BP = bp
    ALT = elev
    TA = Ta[1]
    VEL = vel[1]
    GEOMETRY = shape
    # FLTYPE = 0

    # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
    geom = geometry(o)
    ATOT = geom.AREA
    AV = geom.AV
    AL = geom.AL
    AT = geom.AT

    
    #C     SETTING THE CHARACTERISTIC DIMENSION FOR NUSSELT-REYNOLDS CORRELATIONS
    D = AL

    BETA = 1 / (TA + 273)
    CP = 1.0057E+3
    G = 9.80665

    #C     CONVECTIVE AREA CALCULATION = TOTAL AREA - VENTRAL AREA IN CONTACT WITH SUBSTRATE
    CONVAR = ATOT - AV - AT


    #C     USING ALTITUDE TO COMPUTE BP (SEE DRYAIR LISTING)
    #C      BP=0.0
    DB = TA

    #C     GET THERMAL PROPERTIES OF DRY AIR AT CURRENT TEMP AND PRESSURE
    dair = dryair(db=DB, bp=BP, alt=ALT)
    PATMOS = dair.patmos
    DENSTY = dair.densty
    VISDYN = dair.visdyn
    VISKIN = dair.viskin
    DIFVPR = dair.difvpr
    THCOND = dair.thcond
    HTOVPR = dair.htovpr
    TCOEFF = dair.tcoeff
    GGROUP = dair.ggroup

    #   #C     CHECKING TO SEE IF THE FLUID IS WATER, NOT AIR
    #   if(FLTYPE == 1){
    #     WATERPROP.out = WATERPROP(TA)
    #     CP = WATERPROP.out$CP
    #     DENSTY = WATERPROP.out$DENSTY
    #     THCOND = WATERPROP.out$THCOND
    #     VISDYN = WATERPROP.out$VISDYN
    #   }

    #C     COMPUTING PRANDTL AND SCHMIDT NUMBERS
    PR = CP * VISDYN / THCOND
    #   if(FLTYPE == 0){
    #     #C      AIR
    SC = VISDYN / (DENSTY * DIFVPR)
    #   }else{
    #     #C      WATER; NO MEANING
    #     SC = 1
    #   }

    #C     SKIN/AIR TEMPERATURE DIFFERENCE
    DELTAT = Tskin - TA
    if DELTAT == 0
        DELTAT = 0.01
    end

    #C     COMPUTING GRASHOF NUMBER
    GR = ((DENSTY^2) * BETA * G * (D^3) * DELTAT) / (VISDYN^2)
    #C     CORRECTING IF NEGATIVE DELTAT
    GR = abs(GR)

    #C     AVOIDING DIVIDE BY ZERO IN FREE VS FORCED RAYLEI
    if VEL <= 0
        VEL = 0.0001
    end

    #C     REYNOLDS NUMBER
    RE = DENSTY * VEL * D / VISDYN

    #C     CHOOSING FREE OR FORCED CONVECTION
    #C     SIGNIFICANT FREE CONVECTION IF GR/RE ^ 2 .GE. 20.0
    #C     KREITH (1965) P. 358
    GRRE2 = GR / (RE^2)

    #C     *********************  FREE CONVECTION  ********************
    #   if(GEOMETRY == 0){
    #     RAYLEI = GR * PR
    #     ANUFRE = 0.55 * RAYLEI ^ 0.25
    #   }

    if GEOMETRY == 1 || GEOMETRY == 3 || GEOMETRY == 5
        #C      FREE CONVECTION OF A CYLINDER
        #C      FROM P.334 KREITH (1965): MC ADAM'S 1954 RECOMMENDED COORDINATES
        RAYLEI = GR * PR
        if RAYLEI < 1.0E-05
            ANUFRE = 0.4
        elseif RAYLEI < 0.1
            ANUFRE = 0.976 * RAYLEI^0.0784
        elseif RAYLEI <= 100
            ANUFRE = 1.1173 * RAYLEI^0.1344
        elseif RAYLEI < 10000.0
            ANUFRE = 0.7455 * RAYLEI^0.2167
        elseif RAYLEI < 1.0e9
            ANUFRE = 0.5168 * RAYLEI^0.2501
        elseif RAYLEI < 1.0e12
            ANUFRE = 0.5168 * RAYLEI^0.2501
        end
    end

    #   if((GEOMETRY == 2) | (GEOMETRY == 4)){
    #     #C      SPHERE FREE CONVECTION
    #     #C      FROM P.413 BIRD ET AL (1960) TRANSPORT PHENOMENA)
    #     RAYLEI = (GR ^ (1 /4)) * (PR ^ (1 / 3))
    #     ANUFRE = 2 + 0.60 * RAYLEI
    #     if(RAYLEI < 200){
    #       #GO TO 20
    #     }else{
    #       message(paste0(RAYLEI, '(GR ^ 0.25) * (PR ^ 0.333) IS TOO LARGE FOR CORREL.'))
    #     }
    #   }

    #C     CALCULATING THE FREE CONVECTION HEAT TRANSFER COEFFICIENT, HC  (NU=HC*D/KAIR)
    HCFREE = (ANUFRE * THCOND) / D

    #C     CALCULATING THE SHERWOOD NUMBER FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SHFREE = ANUFRE * (SC / PR)^(1 / 3)
    #C     CALCULATING THE MASS TRANSFER COEFFICIENT FROM THE SHERWOOD NUMBER
    HDFREE = SHFREE * DIFVPR / D
    #C     CALCULATING THE CONVECTIVE HEAT LOSS AT THE SKIN
    QFREE = HCFREE * CONVAR * (Tskin - TA)

    #C     *******************  FORCED CONVECTION  *********************
    #   if(GEOMETRY == 0){
    #     ANU = 0.102 * RE ^ 0.675 * PR ^ (1 / 3)
    #   }

    #   if(GEOMETRY == 1){
    #     #C      FORCED CONVECTION OF A CYLINDER
    #     #C      ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954)
    #     if(RE < 4){
    #       ANU = 0.891 * RE ^ 0.33
    #     }else{
    #     if(RE < 40){
    #       ANU = 0.821 * RE ^ 0.385
    #     }else{
    #     if(RE < 4000.){
    #       ANU = 0.615 * RE ^ 0.466
    #     }else{
    #     if(RE < 40000.){
    #       ANU = 0.174 * RE ^ 0.618
    #     }else{
    #     if(RE < 400000.){
    #       ANU = 0.0239 * RE ^ 0.805
    #     }}}}}
    #   }

    if GEOMETRY > 1
        #C      FORCED CONVECTION IN SPHERE
        #C       ANU=0.34 * RE ^ 0.24 ! ORIGINAL RELATION
        ANU = 0.35 * RE^0.6 # FROM McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
    end

    #C     FORCED CONVECTION FOR ANIMAL

    # if(GEOMETRY == 4){
    #   #C       ***********************FROG******************************
    #   #C       C.R. TRACY'S LEOPARD FROGS - ECOL. MONOG. 1976 V. 46(3)
    #   #C       CHECKING FOR OUT OF BOUNDS VALUES
    #   if(RE < 80){
    #     message(paste0(' RE, ',RE,',TOO SMALL FOR FROG ANCORR'))
    #   }else{
    #     if (RE > 40000){
    #       message(paste0(' RE, ',RE,',TOO LARGE FOR FROG ANCORR'))
    #     }
    #   }
    #   #C       COMPUTING NUSSELT AND SHERWOOD NUMBERS
    #   if(RE <= 2000){
    #     ANU = 0.88 * RE ^ 0.5
    #   }else{
    #     ANU = 0.258 * RE ^ 0.667
    #   }
    # }

    #C     FORCED CONVECTION FOR ANIMAL
    HCFORC = ANU * THCOND / D # HEAT TRANFER COEFFICIENT
    SHFORC = ANU * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HDFORC = SHFORC * DIFVPR / D # MASS TRANSFER COEFFICIENT
    #C     USING BIRD, STEWART, & LIGHTFOOT'S MIXED CONVECTION FORMULA (P. 445, TRANSPORT PHENOMENA, 2002)
    NUTOTAL = (ANUFRE^3 + ANU^3)^(1 / 3)
    HC = NUTOTAL * (THCOND / D) # MIXED CONVECTION HEAT TRANSFER
    QCONV = HC * CONVAR * (Tskin - TA) # TOTAL CONVECTIVE TRANSFER
    #C     CALCULATING THE SHERWOOD NUMBERs FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SH = NUTOTAL * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HD = SH * DIFVPR / D # MASS TRANSFER COEFFICIENT
    return (CONVAR=CONVAR, QCONV=QCONV, HC=HC, HD=HD, SH=SH, QFREE=QFREE, HCFREE=HCFREE, HCFORC=HCFORC, SHFREE=SHFREE, SHFORC=SHFORC, HDFREE=HDFREE, HDFORC=HDFORC)
end



# for zero solver
function convection(Tx, 
    o::fixedParams,
    e::fixedEnvParams,
    v::envVars
)

    Tskin = Tx
    
    @unpack shape = o
    @unpack bp, elev = e
    @unpack Ta, vel = v

    #rename variables: clean the code!
    BP = bp
    ALT = elev
    TA = Ta[1]
    VEL = vel[1]
    GEOMETRY = shape
    # FLTYPE = 0

    # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
    geom = geometry(o)
    ATOT = geom.AREA
    AV = geom.AV
    AL = geom.AL
    AT = geom.AT

    
    #C     SETTING THE CHARACTERISTIC DIMENSION FOR NUSSELT-REYNOLDS CORRELATIONS
    D = AL

    BETA = 1 / (TA + 273)
    CP = 1.0057E+3
    G = 9.80665

    #C     CONVECTIVE AREA CALCULATION = TOTAL AREA - VENTRAL AREA IN CONTACT WITH SUBSTRATE
    CONVAR = ATOT - AV - AT


    #C     USING ALTITUDE TO COMPUTE BP (SEE DRYAIR LISTING)
    #C      BP=0.0
    DB = TA

    #C     GET THERMAL PROPERTIES OF DRY AIR AT CURRENT TEMP AND PRESSURE
    dair = dryair(db=DB, bp=BP, alt=ALT)
    PATMOS = dair.patmos
    DENSTY = dair.densty
    VISDYN = dair.visdyn
    VISKIN = dair.viskin
    DIFVPR = dair.difvpr
    THCOND = dair.thcond
    HTOVPR = dair.htovpr
    TCOEFF = dair.tcoeff
    GGROUP = dair.ggroup

    #   #C     CHECKING TO SEE IF THE FLUID IS WATER, NOT AIR
    #   if(FLTYPE == 1){
    #     WATERPROP.out = WATERPROP(TA)
    #     CP = WATERPROP.out$CP
    #     DENSTY = WATERPROP.out$DENSTY
    #     THCOND = WATERPROP.out$THCOND
    #     VISDYN = WATERPROP.out$VISDYN
    #   }

    #C     COMPUTING PRANDTL AND SCHMIDT NUMBERS
    PR = CP * VISDYN / THCOND
    #   if(FLTYPE == 0){
    #     #C      AIR
    SC = VISDYN / (DENSTY * DIFVPR)
    #   }else{
    #     #C      WATER; NO MEANING
    #     SC = 1
    #   }

    #C     SKIN/AIR TEMPERATURE DIFFERENCE
    DELTAT = Tskin - TA
    if DELTAT == 0
        DELTAT = 0.01
    end

    #C     COMPUTING GRASHOF NUMBER
    GR = ((DENSTY^2) * BETA * G * (D^3) * DELTAT) / (VISDYN^2)
    #C     CORRECTING IF NEGATIVE DELTAT
    GR = abs(GR)

    #C     AVOIDING DIVIDE BY ZERO IN FREE VS FORCED RAYLEI
    if VEL <= 0
        VEL = 0.0001
    end

    #C     REYNOLDS NUMBER
    RE = DENSTY * VEL * D / VISDYN

    #C     CHOOSING FREE OR FORCED CONVECTION
    #C     SIGNIFICANT FREE CONVECTION IF GR/RE ^ 2 .GE. 20.0
    #C     KREITH (1965) P. 358
    GRRE2 = GR / (RE^2)

    #C     *********************  FREE CONVECTION  ********************
    #   if(GEOMETRY == 0){
    #     RAYLEI = GR * PR
    #     ANUFRE = 0.55 * RAYLEI ^ 0.25
    #   }

    if GEOMETRY == 1 || GEOMETRY == 3 || GEOMETRY == 5
        #C      FREE CONVECTION OF A CYLINDER
        #C      FROM P.334 KREITH (1965): MC ADAM'S 1954 RECOMMENDED COORDINATES
        RAYLEI = GR * PR
        if RAYLEI < 1.0E-05
            ANUFRE = 0.4
        elseif RAYLEI < 0.1
            ANUFRE = 0.976 * RAYLEI^0.0784
        elseif RAYLEI <= 100
            ANUFRE = 1.1173 * RAYLEI^0.1344
        elseif RAYLEI < 10000.0
            ANUFRE = 0.7455 * RAYLEI^0.2167
        elseif RAYLEI < 1.0e9
            ANUFRE = 0.5168 * RAYLEI^0.2501
        elseif RAYLEI < 1.0e12
            ANUFRE = 0.5168 * RAYLEI^0.2501
        end
    end

    #   if((GEOMETRY == 2) | (GEOMETRY == 4)){
    #     #C      SPHERE FREE CONVECTION
    #     #C      FROM P.413 BIRD ET AL (1960) TRANSPORT PHENOMENA)
    #     RAYLEI = (GR ^ (1 /4)) * (PR ^ (1 / 3))
    #     ANUFRE = 2 + 0.60 * RAYLEI
    #     if(RAYLEI < 200){
    #       #GO TO 20
    #     }else{
    #       message(paste0(RAYLEI, '(GR ^ 0.25) * (PR ^ 0.333) IS TOO LARGE FOR CORREL.'))
    #     }
    #   }

    #C     CALCULATING THE FREE CONVECTION HEAT TRANSFER COEFFICIENT, HC  (NU=HC*D/KAIR)
    HCFREE = (ANUFRE * THCOND) / D

    #C     CALCULATING THE SHERWOOD NUMBER FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SHFREE = ANUFRE * (SC / PR)^(1 / 3)
    #C     CALCULATING THE MASS TRANSFER COEFFICIENT FROM THE SHERWOOD NUMBER
    HDFREE = SHFREE * DIFVPR / D
    #C     CALCULATING THE CONVECTIVE HEAT LOSS AT THE SKIN
    QFREE = HCFREE * CONVAR * (Tskin - TA)

    #C     *******************  FORCED CONVECTION  *********************
    #   if(GEOMETRY == 0){
    #     ANU = 0.102 * RE ^ 0.675 * PR ^ (1 / 3)
    #   }

    #   if(GEOMETRY == 1){
    #     #C      FORCED CONVECTION OF A CYLINDER
    #     #C      ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954)
    #     if(RE < 4){
    #       ANU = 0.891 * RE ^ 0.33
    #     }else{
    #     if(RE < 40){
    #       ANU = 0.821 * RE ^ 0.385
    #     }else{
    #     if(RE < 4000.){
    #       ANU = 0.615 * RE ^ 0.466
    #     }else{
    #     if(RE < 40000.){
    #       ANU = 0.174 * RE ^ 0.618
    #     }else{
    #     if(RE < 400000.){
    #       ANU = 0.0239 * RE ^ 0.805
    #     }}}}}
    #   }

    if GEOMETRY > 1
        #C      FORCED CONVECTION IN SPHERE
        #C       ANU=0.34 * RE ^ 0.24 ! ORIGINAL RELATION
        ANU = 0.35 * RE^0.6 # FROM McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
    end

    #C     FORCED CONVECTION FOR ANIMAL

    # if(GEOMETRY == 4){
    #   #C       ***********************FROG******************************
    #   #C       C.R. TRACY'S LEOPARD FROGS - ECOL. MONOG. 1976 V. 46(3)
    #   #C       CHECKING FOR OUT OF BOUNDS VALUES
    #   if(RE < 80){
    #     message(paste0(' RE, ',RE,',TOO SMALL FOR FROG ANCORR'))
    #   }else{
    #     if (RE > 40000){
    #       message(paste0(' RE, ',RE,',TOO LARGE FOR FROG ANCORR'))
    #     }
    #   }
    #   #C       COMPUTING NUSSELT AND SHERWOOD NUMBERS
    #   if(RE <= 2000){
    #     ANU = 0.88 * RE ^ 0.5
    #   }else{
    #     ANU = 0.258 * RE ^ 0.667
    #   }
    # }

    #C     FORCED CONVECTION FOR ANIMAL
    HCFORC = ANU * THCOND / D # HEAT TRANFER COEFFICIENT
    SHFORC = ANU * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HDFORC = SHFORC * DIFVPR / D # MASS TRANSFER COEFFICIENT
    #C     USING BIRD, STEWART, & LIGHTFOOT'S MIXED CONVECTION FORMULA (P. 445, TRANSPORT PHENOMENA, 2002)
    NUTOTAL = (ANUFRE^3 + ANU^3)^(1 / 3)
    HC = NUTOTAL * (THCOND / D) # MIXED CONVECTION HEAT TRANSFER
    QCONV = HC * CONVAR * (Tskin - TA) # TOTAL CONVECTIVE TRANSFER
    #C     CALCULATING THE SHERWOOD NUMBERs FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    SH = NUTOTAL * (SC / PR)^(1 / 3) # SHERWOOD NUMBER
    HD = SH * DIFVPR / D # MASS TRANSFER COEFFICIENT
    return (CONVAR=CONVAR, QCONV=QCONV, HC=HC, HD=HD, SH=SH, QFREE=QFREE, HCFREE=HCFREE, HCFORC=HCFORC, SHFREE=SHFREE, SHFORC=SHFORC, HDFREE=HDFREE, HDFORC=HDFORC)
end