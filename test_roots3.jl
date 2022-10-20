
using Parameters, Roots

include("aux_functions.jl")
include("organism.jl")
include("environment.jl")
include("geometry.jl")
include("qmet.jl")
include("qrin.jl")
include("qrout.jl")
include("qsolar.jl")
include("qcond.jl")
include("convection.jl")
include("skinevap.jl")
include("resp.jl")

anim = organismParams(αa = 0.9)
env = envParams()
# vars = envVariables(QSOLR = [400], Tsky=[10], Tsub = [30])
vars = envVariables(Ta=[20])
state = stateVars()

Qrin = qrin(anim, env, vars)
Qsolar = qsolar(anim, env, vars, 800).Qsolar

function energy_bal(Tx,
    o::fixedParams = anim,
    e::fixedEnvParams = env,
    v::envVars = vars,
    s::stateVariables = state,
    Qrin = Qrin,
    Qsolar = Qsolar)
    
    Qmet = qmet(Tx, o)

    # Qout
    Qrout = qrout(Tx, o)
    Qcond = qcond(Tx, o, v)
    Qconv = convection(Tx, o, e, v).QCONV
    Qsevap = skinevap(Tx, o, e, v).QSEVAP
    Qresp = resp(Tx, o, e, v, Qmet).Qresp

    (Qsolar + Qrin + Qmet) - (Qrout + Qcond + Qconv + Qsevap + Qresp)
    # (Qsolar + Qrin + Qmet) - (Qrout + Qcond + Qconv + Qsevap)
    # (Qsolar + Qrin + Qmet) - (Qcond + Qconv + Qsevap)
    # (Qsolar + Qrin + Qmet) - (Qrout + Qcond + Qconv + Qsevap)
    
end


find_zero(energy_bal, (-50, 100), Bisection())
