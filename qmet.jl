


function qmet(;
    Tskin::Real,
    Ww_g::Real, 
    M1::Float64, 
    M2::Float64, 
    M3::Float64
    )

    #Tc = XTRY

    if Tskin > 50 
        Qmet = 0.0056 * 10^(M3 * 50) * M1 * Ww_g^M2
    elseif Tskin >= 1 
        Qmet = 0.0056 * 10^(M3 * Tskin) * M1 * Ww_g^M2
    else 
        Qmet = 0.01
    end

    return Qmet
    
end

function qmet(o::fixedParams, Tskin::Real)

    @unpack M1, M2, M3, Ww_g = o

    #Tc = XTRY

    if Tskin > 50 
        Qmet = 0.0056 * 10^(M3 * 50) * M1 * Ww_g^M2
    elseif Tskin >= 1 
        Qmet = 0.0056 * 10^(M3 * Tskin) * M1 * Ww_g^M2
    else 
        Qmet = 0.01
    end

    return Qmet
    
end


# for the zero solver
function qmet(Tx, o::fixedParams)

    @unpack M1, M2, M3, Ww_g = o

    Tc = Tx

    if Tc > 50 
        Qmet = 0.0056 * 10^(M3 * 50) * M1 * Ww_g^M2
    elseif Tc >= 1 
        Qmet = 0.0056 * 10^(M3 * Tc) * M1 * Ww_g^M2
    else 
        Qmet = 0.01
    end

    return Qmet
    
end

# ##### TESTING 

# Ww_g = 40
# Tskin = 25
# M_1 = 0.013
# M_2 = 0.8
# M_3 = 0.038


# @show qmet(Ww_g = Ww_g, 
# Tskin = Tskin,
#     M1 = M_1, 
#     M2 = M_2, 
#     M3 = M_3)