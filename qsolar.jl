


function qsolar(;
    Atot::Real, 
    Asil::Real, 
    Av::Real, 
    At::Real, 
    αa::Float64, # ABSAN
    αsub::Float64, # ABSSB
    FAsky::Real, 
    FAsub::Real, 
    FAobj::Real, 
    zen::Real, 
    QSOLR::Real, # incoming solar radiation?
    Pdif::Float64, # proportion of solar radiation that is diffuse (fractional, 0-1)
    shade::Real, 
    postur::Int
)

    zenith = zen * 180 / π

    if zenith < 90 
        Qnorm = (QSOLR / cos(zen))
        if postur != 1 
            Qnorm = QSOLR
        end
        if Qnorm > 1367 
            Qnorm = 1367
        end
        if zenith >= 90 
            Qnorm = 0
        end
        Qs_dir = αa * Asil * (1 - Pdif) * Qnorm * (1 - (shade/100)) # direct solar radiation
    else 
        Qs_dir = 0
        Qnorm = 0
    end

    Qs_sky = αa * FAsky * (Atot - At/2) * Pdif * Qnorm * (1 - (shade/100))
    Qs_sub = αa * FAsub * (Atot - Av - At/2) * (1 - αsub) * QSOLR * (1 - (shade/100))
    Qs_obj = αa * FAobj * (Atot - At/2) * Pdif * Qnorm

    Qs_diffuse = Qs_sky + Qs_sub + Qs_obj

    Qsolar = Qs_dir + Qs_diffuse
    
    return (Qsolar = Qsolar, Qs_diffuse = Qs_diffuse, Qs_sub = Qs_sub, Qs_sky = Qs_sky, Qs_obj = Qs_obj)

end





function qsolar(o::fixedParams,
    e::fixedEnvParams,
    v::envVars,
    QSOLR::Real # incoming solar radiation?
)

    @unpack αa, FAsky, FAsub, FAobj, postur = o
    @unpack αsub, Pdif, shade = e
    @unpack zen = v

    geom = geometry(o)
    Atot = geom.AREA
    Av = geom.AV
    At = geom.AT

    # set silhouette area
    if postur == 1
        Asil = geom.ASILN
    elseif postur == 0
        Asil = (geom.ASILN + geom.ASILP) / 2
    end

    zenith = zen[1]
    zen = zen[1] * π / 180

    if zenith < 90 
        Qnorm = (QSOLR / cos(zen))
        if postur != 1 
            Qnorm = QSOLR
        end
        if Qnorm > 1367 
            Qnorm = 1367
        end
        if zenith >= 90 
            Qnorm = 0
        end
        Qs_dir = αa * Asil * (1 - Pdif) * Qnorm * (1 - (shade/100)) # direct solar radiation
    else 
        Qs_dir = 0
        Qnorm = 0
    end

    Qs_sky = αa * FAsky * (Atot - At/2) * Pdif * Qnorm * (1 - (shade/100))
    Qs_sub = αa * FAsub * (Atot - Av - At/2) * (1 - αsub) * QSOLR * (1 - (shade/100))
    Qs_obj = αa * FAobj * (Atot - At/2) * Pdif * Qnorm

    Qs_diffuse = Qs_sky + Qs_sub + Qs_obj

    Qsolar = Qs_dir + Qs_diffuse
    
    return (Qsolar = Qsolar, Qs_diffuse = Qs_diffuse, Qs_sub = Qs_sub, Qs_sky = Qs_sky, Qs_obj = Qs_obj)

end

# #### TESTING

# ATOT = 0.01325006
# ASIL = 0.004718043
# AV = 0.001325006
# AT = 0
# ABSAN = 0.85
# ABSSB = 0.8 
# FATOSK = 0.4
# FATOSB = 0.4
# FATOBJ = 0
# ZEN = 20 * pi/180
# QSOLR = 1000
# PDIF = 0.1
# SHADE = 0
# POSTUR = 1

# @show qsolar(Atot = ATOT, 
#     Asil = ASIL, 
#     Av = AV, 
#     At = AT, 
#     αa = ABSAN, # ABSAN
#     αsub = ABSSB, # ABSSB
#     FAsky = FATOSK, 
#     FAsub = FATOSB, 
#     FAobj = FATOBJ, 
#     zen = ZEN, 
#     QSOLR = QSOLR, # incoming solar radiation?
#     Pdif = PDIF, # proportion of solar radiation that is diffuse (fractional, 0-1)
#     shade = SHADE, 
#     postur = POSTUR)
