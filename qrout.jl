

function qrout(;
    Tskin,
    Atot,
    Av,
    FAsky,
    FAsub,
    εa)

    Ts = c2k(Tskin)
    
    Qr2sky = Atot * FAsky * εa * σ * Ts ^ 4
    Qr2sub = (Atot - Av) * FAsub * εa * σ * Ts ^ 4
    Qr_out = Qr2sky + Qr2sub
    
    return Qr_out
end



function qrout(o::fixedParams,
    s::stateVariables
)

    @unpack FAsky, FAsub, εa = o
    @unpack Tskin = s
    
    # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
    geom = geometry(o)
    Atot = geom.AREA
    Av = geom.AV
    
    Ts = c2k(Tskin)
    
    Qr2sky = Atot * FAsky * εa * σ * Ts ^ 4
    Qr2sub = (Atot - Av) * FAsub * εa * σ * Ts ^ 4
    Qr_out = Qr2sky + Qr2sub
    
    return Qr_out
end


# for zero solver
function qrout(Tx,
    o::fixedParams
)

    Tskin = Tx
    
    @unpack FAsky, FAsub, εa = o
    
    # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
    geom = geometry(o)
    Atot = geom.AREA
    Av = geom.AV
    
    Ts = c2k(Tskin)
    
    Qr2sky = Atot * FAsky * εa * σ * Ts ^ 4
    Qr2sub = (Atot - Av) * FAsub * εa * σ * Ts ^ 4
    Qr_out = Qr2sky + Qr2sub
    
    return Qr_out
end

# TSKIN = 25.1
# ATOT = 0.01325006
# AV = 0.001325006
# AT = 0
# FATOSK = 0.4
# FATOSB = 0.4
# EMISAN = 0.95


# qrout(Tskin = TSKIN,
#     Atot = ATOT,
#     Av = AV,
#     FAsky = FATOSK,
#     FAsub = FATOSB,
#     εa = EMISAN)