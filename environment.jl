using Parameters

abstract type fixedEnvParams end

# parameters that dont change during simulation
@with_kw struct envParams <: fixedEnvParams
    elev::Float64 = 0
    αsub::Float64 = 0.2
    εsub::Float64 = 1
    εsky::Float64 = 1
    bp::Float64 = 101325
    fluid::Int = 0
    O2gas::Float64 = 20.95
    CO2gas::Float64 = 0.03
    N2gas::Float64 = 79.02
    Pdif::Float64 = 0.1
    shade::Float64 = 0
end




abstract type envVars end

@with_kw mutable struct envVariables <: envVars
    QSOLR::Vector{Float64} = [1000]
    zen::Vector{Float64} = [20]
    Ta::Vector{Float64} = [20]
    TGRD::Vector{Float64} = [30]
    Tsub::Vector{Float64} = [30]
    Tsky::Vector{Float64} = [-5]
    vel::Vector{Float64} = [1]
    rh::Vector{Float64} = [5]
    ksub::Vector{Float64} = [0.5]
end

# environ = envVariables(QSOLR=[800])
# environ = envVariables(QSOLR=0:1000)