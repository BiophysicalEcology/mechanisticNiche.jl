


# Trying to figure out the best way of setting parameters

using Parameters

abstract type fixedParams end

# parameters that dont change during simulation
@with_kw struct organismParams <: fixedParams
    Ww_g::Float64 = 40
    αa::Float64 = 0.85
    εa::Float64 = 0.95
    rho_body::Float64 = 1000
    FAsky::Float64 = 0.4
    FAsub::Float64 = 0.4
    FAobj::Float64 = 0.0
    shape::Int = 3
    shape_abc::Vector{Float64} = [1, 3, 2/3]
    custom_shape::Vector{Float64} = [10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743]
    pct_cond::Float64 = 10
    pct_touch::Float64 = 0
    postur::Int = 1
    k_flesh::Float64 = 0.5
    M1::Float64 = 0.013
    M2::Float64 = 0.8
    M3::Float64 = 0.038
    pct_wet::Float64 = 0.1
    pct_eyes::Float64 = 0
    pct_mouth::Float64 = 0
    ψbody::Float64 = -707
    pant::Float64 = 1
    EEf::Float64 = 20 # EXTREF extraction efficiency
    F_O2::Float64 = 20
    RQ::Float64 = 0.8
    ΔT::Float64 = 0.1
    RINSUL::Float64 = 0 # radius of insulation, not used yet
end



abstract type stateVariables end

@with_kw mutable struct stateVars <: stateVariables
    Tc::Float64 = 25
    Tskin::Float64 = 25.1
    Tlung::Float64 = 25
end

# pars = organismParams()
# isa(pars, fixedParams)
# isa(pars, organismParams)
# @unpack Ww_g, alpha, pct_wet = pars;
# Ww_g


# # parameter name translations from input arguments, and value conversions

# AMASS = Ww_g / 1000 # animal wet weight (kg)
# SKINW = pct_wet / 100 # fractional skin wetness
# PEYES = pct_eyes / 100 # fractional of surface area that is wet eyes
# PMOUTH = pct_mouth / 100 # fraction of surface area that is wet mouth
# SKINT = pct_touch / 100 # fraction of surface area that is touching another individual at the same temperature
# PTCOND = pct_cond / 100 # fraction of surface area conducting to the ground
# # PSI_BODY = psi_body #
# ABSAN = alpha # animal solar absorbtivity
# EXTREF = F_O2 # oxygen extraction efficiency (0-1)
# GEOMETRY = shape # animal shape
# FATOSK = fatosk # radiation configuration factor to sky
# FATOSB = fatosb # radiation configuration factor to ground
# FLSHCOND = k_flesh # flesh thermal conductivity (W/mK)
# ANDENS = rho_body # body density (kg/m3)
# PANT = pantmax # panting modifier, >=1 (or 0 if cutting out respiration)
# EMISAN = epsilon # emissivity of animal skin (0-1)
# DELTAR = delta_air # temperature difference between inspired and expired air
# CUSTOMGEOM = custom_shape # parameters for customised geometry

# # unused parameters
# FATOBJ = 0 # configuration factor to nearby object of different temp to sky and ground (e.g. warm rock, fire), not yet used
# RINSUL = 0 # radius of insulation, not used yet
