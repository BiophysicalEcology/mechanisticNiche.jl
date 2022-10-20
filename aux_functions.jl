
# vaporpres AKA VAPPRS

function vaporpres(db::Any)
    t = db .+ 273.16

    if length(db) > 1
        t[t .<= 273.16] = -9.09718 .* (273.16 ./ t[t .<= 273.16] .- 1) .- 
                        3.56654 .* log10.(273.16 ./ t[t .<= 273.16]) .+ 
                        0.876793 .* (1 .- t[t .<= 273.16] ./ 273.16) .+ 
        log10(6.1071)
        t[t .> 273.16] = -7.90298 .* (373.16 ./ t[t .> 273.16] .- 1) .+ 
                        5.02808 .* log10.(373.16 ./ t[t .> 273.16]) .- 
                        1.3816e-07 .* (10 .^ (11.344 .* (1 .- t[t .> 273.16] ./ 373.16)) .- 1) .+ 
                        0.0081328 .* (10 .^ (-3.49149 .* (373.16 ./ t[t .> 273.16] .- 1)) .- 1) .+ 
                        log10(1013.246)
        esat = (10 .^ t) .* 100
    end
    if length(db) == 1
        if t <= 273.16

            t = -9.09718 * (273.16 / t - 1) - 
                3.56654 * log10(273.16 / t) + 
                0.876793 * (1 - t / 273.16) + 
                log10(6.1071)
        else

            t = -7.90298 * (373.16 / t[t > 273.16] - 1) + 
                5.02808 * log10(373.16 / t) - 
                1.3816e-07 * (10 ^ (11.344 * (1 - t / 373.16)) - 1) + 
                0.0081328 * (10 ^ (-3.49149 * (373.16 / t - 1)) - 1) + 
                log10(1013.246)
        end
        esat = (10 ^ t) * 100
    end
    return esat
end

# vaporpres(collect(-10:10))





## wetair function needed in resp

function wetair(; 
    db::Any, 
    wb::Any = db, 
    rh::Any = 0, 
    dp::Any = 999, 
    bp::Any = 101325)
    
    tk = db .+ 273.15
    esat = vaporpres(db)
    if dp < 999                # check this
        e = vaporpres(dp)
        rh = (e ./ esat) .* 100
    else 
        if min(rh) > -1 
            e = esat .* rh ./ 100
        
        else
            wbd = db .- wb
            wbsat = vaporpres(wb)
            dltae = 0.00066 .* (1 .+ 0.00115 .* wb) .* bp .* wbd
            e = wbsat .- dltae
            rh = (e ./ esat) .* 100
        end
    end

    rw = ((0.62197 .* 1.0053 .* e) ./ (bp .- 1.0053 .* e))
    vd = e .* 0.018016 ./ (0.998 .* 8.31434 .* tk)
    tvir = tk .* ((1 .+ rw ./ (18.016 ./ 28.966)) ./ (1 .+ rw))
    tvinc = tvir .- tk
    denair = 0.0034838 .* bp ./ (0.999 .* tvir)
    cp = (1004.84 .+ (rw .* 1846.4)) ./ (1 .+ rw)
    if min(rh) <= 0
        wtrpot = -999
    else
        wtrpot = 461500 .* tk .* log(rh./100)
    end

    return (e = e, esat = esat, vd = vd, rw = rw, tvinc = tvinc, 
            denair = denair, cp = cp, wtrpot = wtrpot, rh = rh)

end
    

# db = 20
# wb = 0
# rh = 50
# dp = 999
# bp = 101325

# db = collect(20:25)

# wetair(db = db, wb = db, rh = rh, dp = dp, bp = bp).e



## dryair function needed in convection

## NEEDS TO BE CHECKED

function dryair(;
    db::Any,
    bp::Any = 0, 
    alt::Any = 0
    )

    tstd = 273.15
    pstd = 101325
    patmos = pstd * ((1 - (0.0065 * alt/288))^(1/0.190284))
    # bp = rep(bp, length(patmos))
    # bp[bp <= 0] = patmos[bp <= 0]
    densty = bp/(287.04 * (db + tstd))
    visnot = 1.8325e-05
    tnot = 296.16
    c = 120
    visdyn = (visnot * ((tnot + c)/(db + tstd + c))) * (((db + 
        tstd)/tnot)^1.5)
    viskin = visdyn/densty
    difvpr = 2.26e-05 * (((db + tstd)/tstd)^1.81) * (1e+05/bp)
    thcond = 0.02425 + (7.038e-05 * db)
    htovpr = 2501200 - 2378.7 * db
    tcoeff = 1/(db + tstd)
    ggroup = 0.0980616 * tcoeff/(viskin * viskin)
    bbemit = 5.670367e-08 * ((db + tstd)^4)
    emtmax = 0.002897/(db + tstd)
    return (patmos = patmos, densty = densty, visdyn = visdyn, 
        viskin = viskin, difvpr = difvpr, thcond = thcond, htovpr = htovpr, 
        tcoeff = tcoeff, ggroup = ggroup, bbemit = bbemit, emtmax = emtmax)

end


