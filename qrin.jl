

# julia version of NicheMapR RADIN function for calculating absorbed long-wave radiation (QIRIN)

const σ = 5.669e-08 # Boltzman's constant

c2k(tempC) = tempC + 273.15


function qrin(; 
    Atot::Real, 
    Av::Real, 
    At::Real, 
    FAsky::Real,
    FAsub::Real, 
    FAobj::Real, 
    εa::Float64, # EMISAN 
    εsub::Float64, # EMISSB 
    εsky::Float64, # EMISSK 
    Tsky::Real, 
    Tsub::Real # TGRD
)

    
    TKsky = c2k(Tsky)
    TKsub = c2k(Tsub)
    TKobj = c2k(Tsub)
    
    Fsky = FAsky - FAobj
    if Fsky < 0 
        Fsky = 0
    end

    Qr_sky = εa * Fsky * (Atot - At / 2) * εsky * σ * TKsky^4
    Qr_sub = εa * FAsub * (Atot - Av - At / 2) * εsub * σ * TKsub^4
    Qr_obj = εa * FAobj * (Atot - At / 2) * εsub * σ * TKobj^4
    
    Qr_in = Qr_sky + Qr_sub + Qr_obj

    return Qr_in
  
end

function qrin(o::fixedParams,
    e::fixedEnvParams,
    v::envVars
)

    @unpack FAsky, FAsub, FAobj, εa = o
    @unpack εsky, εsub = e
    @unpack Tsky, Tsub = v
    
    TKsky = c2k(Tsky[1])
    TKsub = c2k(Tsub[1])
    TKobj = c2k(Tsub[1])
    
    Fsky = FAsky - FAobj
    if Fsky < 0 
        Fsky = 0
    end

    # or maybe in a later version include a geom object as input (only call geom ~once per time-step)
    geom = geometry(o)
    Atot = geom.AREA
    Av = geom.AV
    At = geom.AT

    Qr_sky = εa * Fsky * (Atot - At / 2) * εsky * σ * TKsky^4
    Qr_sub = εa * FAsub * (Atot - Av - At / 2) * εsub * σ * TKsub^4
    Qr_obj = εa * FAobj * (Atot - At / 2) * εsub * σ * TKobj^4
    
    Qr_in = Qr_sky + Qr_sub + Qr_obj

    return Qr_in
  
end


