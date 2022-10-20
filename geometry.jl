



function geometry(; AMASS = 0.04,
    GEOMETRY = 3,
    SHP = [1, 3, 2 / 3],
    CUSTOMGEOM = [10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743],
    ANDENS = 1000,
    SKINW = 0.001,
    SKINT = 0,
    RINSUL = 0,
    PTCOND = 0.1,
    PMOUTH = 0.05,
    PANT = 1)

    VOL = AMASS / ANDENS
    R = ((3 / 4) * VOL / π) ^ (1 / 3)
    R1 = R - RINSUL
    D = 2 * R1
    AL = VOL ^ (1 / 3) #JOHN MITCHELL'S CHARACTERISTIC DIMENSION FOR CONVECTION (1976)

    GMASS =  AMASS * 1000
    ASEMAJR = 0
    BSEMINR = 0
    CSEMINR = 0
    AT = 0

    # #C     FLAT PLATE
    # if(GEOMETRY == 0){
    #     #C      ASSUME A SQUARE FOR THE MOMENT
    #     AHEIT = (VOL / (SHP[2]*SHP[3])) ^ (1 / 3)
    #     AWIDTH = AHEIT*SHP[2]
    #     ALENTH =  AHEIT*SHP[3]
    #     ATOT = ALENTH * AWIDTH * 2 + ALENTH * AHEIT * 2 + AWIDTH * AHEIT * 2
    #     AREA = ATOT
    #     ASILN = ALENTH * AWIDTH
    #     ASILP = AWIDTH * AHEIT
    #     AL = AHEIT
    #     if(AWIDTH <= ALENTH){
    #         AL = AWIDTH
    #     }else{
    #         AL = ALENTH
    #     }
    #     R = ALENTH / 2
    #     AV = ATOT  * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    # }

    # #C     CYLINDER
    # if(GEOMETRY == 1){
    #     R1 = (VOL / (PI * SHP[2] * 2)) ^ (1 / 3)
    #     ALENTH = R1 * SHP[2]
    #     AREA = 2 * PI * R1 ^ 2 + 2 * PI * R1 * ALENTH
    #     VOL = AMASS / ANDENS
    #     AWIDTH = 2 * R1
    #     ASILN = AWIDTH  *  ALENTH
    #     ASILP = PI * R1 ^ 2
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     ATOT = AREA
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    #     AL = ALENTH
    #     TOTLEN = ALENTH
    # }

    # #C     ELLIPSOID
    # if(GEOMETRY == 2){
    #     A = ((3 / 4) * VOL / (PI * SHP[2] * SHP[3])) ^ (1 / 3)
    #     B = A * SHP[2]
    #     C = A * SHP[3]
    #     P = 1.6075
    #     AREA = (4 * PI * (((A ^ P * B ^ P + A ^ P * C ^ P + B ^ P * C ^ P)) / 3) ^ (1 / P))
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     ATOT = AREA
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    #     ASILN = PI * A * C
    #     ASILP = PI * B * C
    #     if(ASILN < ASILP){
    #         ASILN = PI * B * C
    #         ASILP = PI * A * C
    #     }
    #     ASEMAJR = A
    #     BSEMINR = B
    #     CSEMINR = C
    # }

    #C     LIZARD
    if GEOMETRY == 3
        ATOTAL = (10.4713 * GMASS ^ 0.688) / 10000.
        AV = (0.425 * GMASS ^ 0.85) / 10000.
        ATOT = ATOTAL
        VOL = AMASS / ANDENS
        #C      CONDUCTION - RADIATION, ETC DIMENSION. ASSUME L=2D=4R1;
        #C      THEN VOL=PI * R1 ^ 2 * L = 4 * PI * R1 ^ 3
        R1 = (VOL / (4 * π)) ^ (1 / 3)
        #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
        #C      MAX. SILHOUETTE AREA (NORMAL TO THE SUN)
        ASILN = (3.798 * GMASS ^ 0.683) / 10000
        #C      MIN. SILHOUETTE AREA (POINTING TOWARD THE SUN)
        ASILP = (0.694 * GMASS ^ 0.743) / 10000
        AREA = ATOT
        AV = AREA * PTCOND
        AT = AREA * SKINT
        if PANT > 1
            AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
        else
            AEFF = SKINW * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
        end
    end

    # if(GEOMETRY == 4){
    #     #C      AREA OF LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
    #     ATOTAL = (12.79 * GMASS ^ 0.606) / 10000.
    #     AV = (0.425 * GMASS ^ 0.85) / 10000.
    #     ATOT = ATOTAL
    #     #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: EQ'N 11 TRACY 1976
    #     ZEN = 0
    #     PCTN = 1.38171E-06 * ZEN ^ 4-1.93335E-04 * ZEN ^ 3 +
    #         4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    #     ASILN = PCTN * ATOT / 100
    #     ZEN = 90
    #     PCTP = 1.38171E-06 * ZEN ^ 4-1.93335E-04 * ZEN ^ 3 +
    #         4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    #     ASILP = PCTP * ATOT / 100
    #     AREA = ATOT
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND-SKINW * AREA * SKINT
    #     }else{
    #         AEFF = SKINW * AREA-SKINW * AREA * PTCOND-SKINW * AREA * SKINT
    #     }
    # }

    # #C     USER-DEFINED GEOMETRY
    # if(GEOMETRY == 5){
    #     ATOTAL = (CUSTOMGEOM[1] * GMASS ^ CUSTOMGEOM[2]) / 10000
    #     AV = (CUSTOMGEOM[3] * GMASS ^ CUSTOMGEOM[4]) / 10000
    #     ATOT = ATOTAL
    #     VOL = AMASS / ANDENS
    #     #C      CONDUCTION - RADIATION, ETC DIMENSION. ASSUME L = 2D = 4R1;
    #     #C      THEN VOL = PI * R1 ^ 2 * L  =  4 * PI * R1 ^ 3
    #     R1 = (VOL / (4 * PI)) ^ (1 / 3)
    #     #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    #     #C      USER MUST DEFINE MAX. SILHOUETTE AREA (NORMAL TO THE SUN)
    #     ASILN = (CUSTOMGEOM[5] * GMASS ^ CUSTOMGEOM[6]) / 10000
    #     #C      USER MUST DEFINE MIN. SILHOUETTE AREA (POINTING TOWARD THE SUN)
    #     ASILP = (CUSTOMGEOM[7] * GMASS ^ CUSTOMGEOM[8]) / 10000
    #     AREA = ATOT
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
    #     }else{
    #         AEFF = SKINW * AREA - SKINW * AREA * PTCOND - SKINW * AREA * SKINT
    #     }
    # }
  
    return (AREA = AREA, AV = AV, AT = AT, AL = AL, VOL = VOL, R = R, R1 = R1, ASEMAJR = ASEMAJR, BSEMINR = BSEMINR, CSEMINR = CSEMINR, ASILN = ASILN, ASILP = ASILP, AEFF = AEFF)

end



# second behavior of the function

function geometry(o::fixedParams)

    @unpack Ww_g, shape, shape_abc, custom_shape, rho_body, pct_wet = o
    @unpack pct_touch, RINSUL, pct_cond, pct_mouth, pant = o

    AMASS = Ww_g / 1000
    GEOMETRY = shape
    SHP = shape_abc
    CUSTOMGEOM = custom_shape
    ANDENS = rho_body
    SKINW = pct_wet / 100
    SKINT = pct_touch / 100
    PTCOND = pct_cond / 100
    PMOUTH = pct_mouth / 100
    PANT = pant

    VOL = AMASS / ANDENS
    R = ((3 / 4) * VOL / π) ^ (1 / 3)
    R1 = R - RINSUL
    D = 2 * R1
    AL = VOL ^ (1 / 3) #JOHN MITCHELL'S CHARACTERISTIC DIMENSION FOR CONVECTION (1976)

    GMASS =  AMASS * 1000
    ASEMAJR = 0
    BSEMINR = 0
    CSEMINR = 0
    AT = 0

    # #C     FLAT PLATE
    # if(GEOMETRY == 0){
    #     #C      ASSUME A SQUARE FOR THE MOMENT
    #     AHEIT = (VOL / (SHP[2]*SHP[3])) ^ (1 / 3)
    #     AWIDTH = AHEIT*SHP[2]
    #     ALENTH =  AHEIT*SHP[3]
    #     ATOT = ALENTH * AWIDTH * 2 + ALENTH * AHEIT * 2 + AWIDTH * AHEIT * 2
    #     AREA = ATOT
    #     ASILN = ALENTH * AWIDTH
    #     ASILP = AWIDTH * AHEIT
    #     AL = AHEIT
    #     if(AWIDTH <= ALENTH){
    #         AL = AWIDTH
    #     }else{
    #         AL = ALENTH
    #     }
    #     R = ALENTH / 2
    #     AV = ATOT  * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    # }

    # #C     CYLINDER
    # if(GEOMETRY == 1){
    #     R1 = (VOL / (PI * SHP[2] * 2)) ^ (1 / 3)
    #     ALENTH = R1 * SHP[2]
    #     AREA = 2 * PI * R1 ^ 2 + 2 * PI * R1 * ALENTH
    #     VOL = AMASS / ANDENS
    #     AWIDTH = 2 * R1
    #     ASILN = AWIDTH  *  ALENTH
    #     ASILP = PI * R1 ^ 2
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     ATOT = AREA
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    #     AL = ALENTH
    #     TOTLEN = ALENTH
    # }

    # #C     ELLIPSOID
    # if(GEOMETRY == 2){
    #     A = ((3 / 4) * VOL / (PI * SHP[2] * SHP[3])) ^ (1 / 3)
    #     B = A * SHP[2]
    #     C = A * SHP[3]
    #     P = 1.6075
    #     AREA = (4 * PI * (((A ^ P * B ^ P + A ^ P * C ^ P + B ^ P * C ^ P)) / 3) ^ (1 / P))
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     ATOT = AREA
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * (AREA - AV)
    #     }else{
    #         AEFF = SKINW * (AREA - AV)
    #     }
    #     ASILN = PI * A * C
    #     ASILP = PI * B * C
    #     if(ASILN < ASILP){
    #         ASILN = PI * B * C
    #         ASILP = PI * A * C
    #     }
    #     ASEMAJR = A
    #     BSEMINR = B
    #     CSEMINR = C
    # }

    #C     LIZARD
    if GEOMETRY == 3
        ATOTAL = (10.4713 * GMASS ^ 0.688) / 10000.
        AV = (0.425 * GMASS ^ 0.85) / 10000.
        ATOT = ATOTAL
        VOL = AMASS / ANDENS
        #C      CONDUCTION - RADIATION, ETC DIMENSION. ASSUME L=2D=4R1;
        #C      THEN VOL=PI * R1 ^ 2 * L = 4 * PI * R1 ^ 3
        R1 = (VOL / (4 * π)) ^ (1 / 3)
        #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
        #C      MAX. SILHOUETTE AREA (NORMAL TO THE SUN)
        ASILN = (3.798 * GMASS ^ 0.683) / 10000
        #C      MIN. SILHOUETTE AREA (POINTING TOWARD THE SUN)
        ASILP = (0.694 * GMASS ^ 0.743) / 10000
        AREA = ATOT
        AV = AREA * PTCOND
        AT = AREA * SKINT
        if PANT > 1
            AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
        else
            AEFF = SKINW * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
        end
    end

    # if(GEOMETRY == 4){
    #     #C      AREA OF LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
    #     ATOTAL = (12.79 * GMASS ^ 0.606) / 10000.
    #     AV = (0.425 * GMASS ^ 0.85) / 10000.
    #     ATOT = ATOTAL
    #     #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: EQ'N 11 TRACY 1976
    #     ZEN = 0
    #     PCTN = 1.38171E-06 * ZEN ^ 4-1.93335E-04 * ZEN ^ 3 +
    #         4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    #     ASILN = PCTN * ATOT / 100
    #     ZEN = 90
    #     PCTP = 1.38171E-06 * ZEN ^ 4-1.93335E-04 * ZEN ^ 3 +
    #         4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    #     ASILP = PCTP * ATOT / 100
    #     AREA = ATOT
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND-SKINW * AREA * SKINT
    #     }else{
    #         AEFF = SKINW * AREA-SKINW * AREA * PTCOND-SKINW * AREA * SKINT
    #     }
    # }

    # #C     USER-DEFINED GEOMETRY
    # if(GEOMETRY == 5){
    #     ATOTAL = (CUSTOMGEOM[1] * GMASS ^ CUSTOMGEOM[2]) / 10000
    #     AV = (CUSTOMGEOM[3] * GMASS ^ CUSTOMGEOM[4]) / 10000
    #     ATOT = ATOTAL
    #     VOL = AMASS / ANDENS
    #     #C      CONDUCTION - RADIATION, ETC DIMENSION. ASSUME L = 2D = 4R1;
    #     #C      THEN VOL = PI * R1 ^ 2 * L  =  4 * PI * R1 ^ 3
    #     R1 = (VOL / (4 * PI)) ^ (1 / 3)
    #     #C      NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    #     #C      USER MUST DEFINE MAX. SILHOUETTE AREA (NORMAL TO THE SUN)
    #     ASILN = (CUSTOMGEOM[5] * GMASS ^ CUSTOMGEOM[6]) / 10000
    #     #C      USER MUST DEFINE MIN. SILHOUETTE AREA (POINTING TOWARD THE SUN)
    #     ASILP = (CUSTOMGEOM[7] * GMASS ^ CUSTOMGEOM[8]) / 10000
    #     AREA = ATOT
    #     AV = AREA * PTCOND
    #     AT = AREA * SKINT
    #     if(PANT > 1){
    #         AEFF = (SKINW + PMOUTH) * AREA-SKINW * AREA * PTCOND - SKINW * AREA * SKINT
    #     }else{
    #         AEFF = SKINW * AREA - SKINW * AREA * PTCOND - SKINW * AREA * SKINT
    #     }
    # }
  
    return (AREA = AREA, AV = AV, AT = AT, AL = AL, VOL = VOL, R = R, R1 = R1, ASEMAJR = ASEMAJR, BSEMINR = BSEMINR, CSEMINR = CSEMINR, ASILN = ASILN, ASILP = ASILP, AEFF = AEFF)

end

# @show geometry()