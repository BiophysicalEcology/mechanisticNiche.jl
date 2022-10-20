

function qcond(;
    Av = 0.001325006,
    dep = 0.025,
    Tskin = 25,
    Tsub = 10,
    ksub = 0.1 # previously SUBTK. Thermal conductivity (k) of the substrate (sub)?
)

    Qcond = Av * (ksub / dep) * (Tskin - Tsub) 

    return Qcond
end



function qcond(o::fixedParams,
    v::envVars,
    s::stateVariables,
    dep::Real = 0.025
)

    @unpack Tsub, ksub = v
    @unpack Tskin = s

    Tsub = Tsub[1]
    ksub = ksub[1]

    geom = geometry(o)
    Av = geom.AV

    Qcond = Av * (ksub / dep) * (Tskin - Tsub) 

    return Qcond
end


# for zero solver
function qcond(Tx,
    o::fixedParams,
    v::envVars,
    dep::Real = 0.025
)

    Tskin = Tx
    @unpack Tsub, ksub = v
    
    Tsub = Tsub[1]
    ksub = ksub[1]

    geom = geometry(o)
    Av = geom.AV

    Qcond = Av * (ksub / dep) * (Tskin - Tsub) 

    return Qcond
end


# AV = 0.001325006
# DEP = 0.025
# TSKIN = 25
# TSUBST = 10
# SUBTK = 0.1

# qcond(Av = AV,
# dep = DEP,
# Tskin = TSKIN,
# Tsub = TSUBST,
# SUBTK = SUBTK)
